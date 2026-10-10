using System.Reflection;
using System.Threading;

namespace PatchedConicFixes
{
    /// <summary>
    ///     Replacement for the stock Orbit.FindClosestPoints delegate.  The stock implementation (Orbit._FindClosestPoints) runs
    ///     the recursive sampling search in the Targeting class, which can miss minima and can fail to terminate.  This routes every
    ///     caller of that API, including any other mods, through Baluev's solver instead.
    ///     The API can only return two minima, but two orbits can have up to four (and FindAllMinima() can return up to six).  Callers
    ///     of this API pass in no time information, so the "best" two are the two with the smallest geometrical distance, which is the
    ///     global MOID and the next most likely encounter site.  The stock callers that can make a better choice using time information
    ///     are patched to use InterceptDisplay.SelectDisplayIntercepts() instead.
    /// </summary>
    public static class FindClosestPointsAdapter
    {
        private static int _fallbackLogged;

        /// <summary>
        ///     Points Orbit.FindClosestPoints at FindClosestPoints() below, unless another mod has already replaced it.
        /// </summary>
        public static void Install()
        {
            MethodInfo current = Orbit.FindClosestPoints.Method;

            if (current.DeclaringType == typeof(FindClosestPointsAdapter))
                return;

            if (current != typeof(Orbit).GetMethod(nameof(Orbit._FindClosestPoints)))
            {
                Logger.Print($"Orbit.FindClosestPoints was already replaced by {current.DeclaringType}.{current.Name}, leaving it alone");
                return;
            }

            Orbit.FindClosestPoints = FindClosestPoints;
            Logger.Print("Orbit.FindClosestPoints now uses the Baluev MOID solver");
        }

        /// <summary>
        ///     Drop-in replacement for Orbit._FindClosestPoints().  The outputs match the stock API:
        ///     CD/FFp/FFs and CCD/SFp/SFs are the distance, the true anomaly on p and the true anomaly on s of the
        ///     first and second minima.  If only one minimum is found it is copied into both slots, and if none are found the
        ///     outputs are left untouched.  The epsilon, maxIterations and iterationCount parameters are unused (as they are in stock).
        /// </summary>
        /// <returns>The number of minima found (0, 1 or 2)</returns>
        public static unsafe int FindClosestPoints(Orbit p, Orbit s, ref double CD, ref double CCD, ref double FFp, ref double FFs, ref double SFp, ref double SFs, double epsilon, int maxIterations, ref int iterationCount)
        {
            if (GameSettings.LEGACY_ORBIT_TARGETING)
                return Orbit._FindClosestPoints(p, s, ref CD, ref CCD, ref FFp, ref FFs, ref SFp, ref SFs, epsilon, maxIterations, ref iterationCount);

            Baluev.MoidInfo* info = stackalloc Baluev.MoidInfo[6];

            int num = HarmonyPatches.FindMoidMinima(p, s, info);

            // FindAllMinima() can't solve every case (e.g. exactly parabolic orbits), so let the stock solver have a go at those.
            if (num < 1)
            {
                if (Interlocked.Exchange(ref _fallbackLogged, 1) == 0)
                    Logger.Print($"FindClosestPoints: no minima from the Baluev solver, falling back to the stock solver (p: {HarmonyPatches.OrbitString(p)}, s: {HarmonyPatches.OrbitString(s)})");

                return Orbit._FindClosestPoints(p, s, ref CD, ref CCD, ref FFp, ref FFs, ref SFp, ref SFs, epsilon, maxIterations, ref iterationCount);
            }

            int count = SameSideMinima(p, s, info, num);

            int first = -1, second = -1;

            for (int i = 0; i < count; i++)
            {
                if (first < 0 || info[i].dst < info[first].dst)
                {
                    second = first;
                    first  = i;
                }
                else if (second < 0 || info[i].dst < info[second].dst)
                {
                    second = i;
                }
            }

            return CopyMinima(info, first, second, ref CD, ref CCD, ref FFp, ref FFs, ref SFp, ref SFs);
        }

        /// <summary>
        ///     Removes minima where the two points are on opposite sides of the central body, as the stock solver does.
        /// </summary>
        /// <param name="info">Minima, with ta1 and ta2 resolved.  Compacted in place.</param>
        /// <param name="num">The number of minima in info</param>
        /// <returns>The number of minima left in info</returns>
        public static unsafe int SameSideMinima(Orbit p, Orbit s, Baluev.MoidInfo* info, int num)
        {
            int count = 0;
            for (int i = 0; i < num; i++)
            {
                if (Vector3d.Dot(p.getRelativePositionFromTrueAnomaly(info[i].ta1), s.getRelativePositionFromTrueAnomaly(info[i].ta2)) <= 0)
                    continue;

                info[count++] = info[i];
            }

            return count;
        }

        /// <summary>
        ///     Fills in the outputs of the FindClosestPoints API from the selected minima, matching the stock behavior when fewer than two
        ///     were selected.
        /// </summary>
        /// <param name="info">Minima, with ta1 and ta2 resolved</param>
        /// <param name="first">Index of the first minimum, or -1 if there is none</param>
        /// <param name="second">Index of the second minimum, or -1 if there is none</param>
        /// <returns>The number of minima returned (0, 1 or 2)</returns>
        public static unsafe int CopyMinima(Baluev.MoidInfo* info, int first, int second, ref double CD, ref double CCD, ref double FFp, ref double FFs, ref double SFp, ref double SFs)
        {
            if (first < 0)
                return 0;

            CD  = info[first].dst;
            FFp = info[first].ta1;
            FFs = info[first].ta2;

            if (second < 0)
            {
                CCD = CD;
                SFp = FFp;
                SFs = FFs;
                return 1;
            }

            CCD = info[second].dst;
            SFp = info[second].ta1;
            SFs = info[second].ta2;
            return 2;
        }
    }
}
