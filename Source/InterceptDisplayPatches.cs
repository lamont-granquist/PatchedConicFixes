using System.Collections.Generic;
using System.Reflection;
using System.Reflection.Emit;
using HarmonyLib;

// ReSharper disable ArrangeTypeModifiers
// ReSharper disable UnusedMember.Local
// ReSharper disable UnusedType.Global
// ReSharper disable ArrangeTypeMemberModifiers

namespace PatchedConicFixes
{
    // The Intercept tab of the maneuver node editor.  Redirects its Orbit.FindClosestPoints call (see issue #2).
    [HarmonyPatch(typeof(OrbitUtil), nameof(OrbitUtil.CalculateIntersections))]
    class OrbitUtil_CalculateIntersections
    {
        static IEnumerable<CodeInstruction> Transpiler(IEnumerable<CodeInstruction> instructions) =>
            InterceptDisplay.RedirectFindClosestPoints(instructions, "OrbitUtil.CalculateIntersections");
    }

    // The intersection markers in map view when targeting a vessel.  Redirects its Orbit.FindClosestPoints call (see issue #3).
    [HarmonyPatch(typeof(OrbitTargeter), "UpdateNodesAndVectors")]
    class OrbitTargeter_UpdateNodesAndVectors
    {
        static IEnumerable<CodeInstruction> Transpiler(IEnumerable<CodeInstruction> instructions) =>
            InterceptDisplay.RedirectFindClosestPoints(instructions, "OrbitTargeter.UpdateNodesAndVectors");
    }

    public static class InterceptDisplay
    {
        /// <summary>
        ///     Replaces the single Orbit.FindClosestPoints call in a stock method with a call to SelectDisplayIntercepts().
        ///     If the method doesn't look as expected it is left unmodified (and will then use FindClosestPointsAdapter through the
        ///     delegate).  This deliberately doesn't throw, since that would abort Harmony's PatchAll() and leave other patches unapplied.
        /// </summary>
        public static IEnumerable<CodeInstruction> RedirectFindClosestPoints(IEnumerable<CodeInstruction> instructions, string caller)
        {
            var code = new List<CodeInstruction>(instructions);

            FieldInfo  field  = AccessTools.Field(typeof(Orbit), nameof(Orbit.FindClosestPoints));
            MethodInfo invoke = AccessTools.Method(typeof(Orbit.FindClosestPointsDelegate), nameof(Orbit.FindClosestPointsDelegate.Invoke));
            MethodInfo select = AccessTools.Method(typeof(InterceptDisplay), nameof(SelectDisplayIntercepts));

            int load = -1, call = -1, loads = 0, calls = 0;

            for (int i = 0; i < code.Count; i++)
            {
                if (code[i].LoadsField(field))
                {
                    load = i;
                    loads++;
                }
                else if (code[i].Calls(invoke))
                {
                    call = i;
                    calls++;
                }
            }

            if (loads != 1 || calls != 1 || load > call)
            {
                Logger.Print($"{caller}: expected one call to Orbit.FindClosestPoints but found {loads} loads and {calls} calls, leaving it unpatched");
                return code;
            }

            // The delegate instance is pushed before the arguments, so stop loading it, then call our static method (which has the same
            // parameters) instead of the delegate.  The load becomes a nop rather than being removed so that any labels on it survive.
            code[load].opcode  = OpCodes.Nop;
            code[load].operand = null;
            code[call].opcode  = OpCodes.Call;
            code[call].operand = select;

            return code;
        }

        /// <summary>
        ///     FindClosestPoints() replacement for the UI code that displays the intercepts with a target vessel.  That code can
        ///     only display two intercepts, so of all the minima that Baluev's solver finds this returns the two where the target
        ///     will actually be closest when p passes the minimum, which is what the UI displays.  The UT computed here is the
        ///     same one that the stock callers compute from the returned true anomalies.  Minima outside of p's patch bounds, or with the
        ///     two points on opposite sides of the central body, are skipped, and orbits around different bodies have no intercepts.
        ///     Otherwise the outputs match FindClosestPointsAdapter.
        /// </summary>
        /// <returns>The number of intercepts found (0, 1 or 2)</returns>
        public static unsafe int SelectDisplayIntercepts(Orbit p, Orbit s, ref double CD, ref double CCD, ref double FFp, ref double FFs, ref double SFp, ref double SFs, double epsilon, int maxIterations, ref int iterationCount)
        {
            if (GameSettings.LEGACY_ORBIT_TARGETING)
                return Orbit.FindClosestPoints(p, s, ref CD, ref CCD, ref FFp, ref FFs, ref SFp, ref SFs, epsilon, maxIterations, ref iterationCount);

            // The maneuver node editor walks through patches in other SOIs, where the geometry is meaningless.
            if (!ReferenceEquals(p.referenceBody, s.referenceBody))
                return 0;

            Baluev.MoidInfo* info = stackalloc Baluev.MoidInfo[6];
            double*          sep  = stackalloc double[6];

            int num = HarmonyPatches.FindMoidMinima(p, s, info);

            // Let FindClosestPointsAdapter fall back to the stock solver for anything Baluev can't solve.
            if (num < 1)
                return Orbit.FindClosestPoints(p, s, ref CD, ref CCD, ref FFp, ref FFs, ref SFp, ref SFs, epsilon, maxIterations, ref iterationCount);

            num = FindClosestPointsAdapter.SameSideMinima(p, s, info, num);

            int count = 0;
            for (int i = 0; i < num; i++)
            {
                if (!PatchedConics.TAIsWithinPatchBounds(info[i].ta1, p))
                    continue;

                double ut = p.StartUT + p.GetDTforTrueAnomaly(info[i].ta1, 0.0);

                if (!Statics.IsFinite(ut))
                    continue;

                sep[count]    = (p.getRelativePositionAtUT(ut) - s.getRelativePositionAtUT(ut)).magnitude;
                info[count++] = info[i];
            }

            int first = -1, second = -1;

            for (int i = 0; i < count; i++)
            {
                if (first < 0 || sep[i] < sep[first])
                {
                    second = first;
                    first  = i;
                }
                else if (second < 0 || sep[i] < sep[second])
                {
                    second = i;
                }
            }

            return FindClosestPointsAdapter.CopyMinima(info, first, second, ref CD, ref CCD, ref FFp, ref FFs, ref SFp, ref SFs);
        }
    }
}
