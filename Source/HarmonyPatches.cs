using System;
using System.Threading;
using HarmonyLib;
using UnityEngine;
using static PatchedConicFixes.Statics;

// ReSharper disable ArrangeTypeModifiers
// ReSharper disable UnusedMember.Local
// ReSharper disable UnusedType.Global
// ReSharper disable ArrangeTypeMemberModifiers
// ReSharper disable RedundantAssignment

namespace PatchedConicFixes
{
    // Replaces PatchedConics._CheckEncounter entirely.
    [HarmonyPatch(typeof(PatchedConics), nameof(PatchedConics._CheckEncounter))]
    class PatchedConics__CheckEncounter
    {
        static bool Prefix(Orbit p, Orbit nextPatch, double startEpoch, OrbitDriver sec, CelestialBody targetBody, PatchedConics.SolverParameters pars, bool logErrors, ref bool __result)
        {
            //__result = HarmonyPatches.CheckEncounter(p, nextPatch, startEpoch, sec, targetBody, pars, logErrors);
            __result = HarmonyPatches.CheckEncounterOldMOID(p, nextPatch, startEpoch, sec, targetBody, pars, logErrors);

            return false;
        }
    }

    // Replaces PatchedConics._GetClosestApproach entirely.
    [HarmonyPatch(typeof(PatchedConics), nameof(PatchedConics._GetClosestApproach))]
    class PatchedConics__GetClosestApproach
    {
        static bool Prefix(Orbit p, Orbit s, double startEpoch, double dT, PatchedConics.SolverParameters pars, ref double __result)
        {
            __result = HarmonyPatches.GetClosestApproach(p, s, startEpoch, dT, pars);

            return false;
        }
    }

    // Replaces Orbit._SolveClosestApproach entirely.
    [HarmonyPatch(typeof(Orbit), nameof(Orbit._SolveClosestApproach))]
    class Orbit__SolveClosestApproach
    {
        static bool Prefix(Orbit p, Orbit s, ref double UT, double dT, double threshold, double MinUT, double MaxUT, double epsilon, int maxIterations, ref int iterationCount, ref double __result)
        {
            __result = HarmonyPatches.SolveClosestApproach(p, s, ref UT, dT, threshold, MinUT, MaxUT, epsilon, maxIterations, ref iterationCount);

            return false;
        }
    }

    // Replaces Orbit._SolveSOI entirely.
    [HarmonyPatch(typeof(Orbit), nameof(Orbit._SolveSOI))]
    class Orbit__SolveSOI
    {
        static bool Prefix(Orbit p, Orbit s, ref double UT, double dT, double Rsoi, double MinUT, double MaxUT, double epsilon, int maxIterations, ref int iterationCount, ref bool __result)
        {
            __result = HarmonyPatches.SolveSOI(p, s, ref UT, dT, Rsoi, MinUT, MaxUT, epsilon, maxIterations, ref iterationCount);

            return false;
        }
    }

    public static class HarmonyPatches
    {
        public static string OrbitString(Orbit o) => $"{o.inclination}, {o.eccentricity}, {o.semiMajorAxis}, {o.LAN}, {o.argumentOfPeriapsis}, {o.meanAnomalyAtEpoch}, {o.epoch}";

        /// <summary>
        ///     Determines whether the primary orbit encounters (enters the SOI of) the given celestial body.
        ///     This is the main encounter detection entry point, implementing a multi-stage filter pipeline:
        ///     Stage 1: Pe/Ap prefilter: Checks whether the altitude bands of the two orbits overlap
        ///     (with a buffer for the body's SOI radius). This is omitted when the target is selected.
        ///     Stage 2: Geometric MOID solver: Calls FindClosestPoints to find the one or two points on
        ///     the orbits where the inter-orbit distance function has critical points (local minima). These
        ///     are purely geometric--they identify WHERE on the orbits a close approach could occur,
        ///     independent of when the bodies actually reach those points.
        ///     Stage 3: Validation of the critical points against the patch boundaries.
        ///     Stage 4: Time-domain closest approach: For each valid geometric intercept (up to two),
        ///     runs the time-domain Halley solver (GetClosestApproach) seeded at that intercept to find the
        ///     actual minimum separation accounting for both bodies' motion. If the closest approach is
        ///     within the body's SOI, calls EncountersBody to compute the SOI crossing time and generate
        ///     the next orbit patch.
        /// </summary>
        /// <param name="p">
        ///     The primary orbit (vessel patch being tested). Modified in place with encounter
        ///     data if an encounter is found (UTappr, ClAppr, UTsoi, EndUT, patchEndTransition, etc).
        /// </param>
        /// <param name="nextPatch">
        ///     The orbit object to populate with the post-encounter trajectory if an
        ///     SOI transition occurs.  Caller must provide a new Orbit.
        /// </param>
        /// <param name="startEpoch">The UT at the start of this patch.</param>
        /// <param name="sec">The OrbitDriver of the celestial body being tested for encounter.</param>
        /// <param name="targetBody">
        ///     The player's currently targeted body, or null. Used to compute
        ///     approach markers even when no SOI encounter occurs.
        /// </param>
        /// <param name="pars">Solver parameters (iteration limits, epsilon, etc).</param>
        /// <param name="logErrors">If false, suppresses debug logging. Used by background thread callers.</param>
        /// <returns>
        ///     True if the primary orbit enters the body's SOI (and nextPatch has been populated
        ///     with the post-transition orbit), false otherwise.
        /// </returns>
        public static unsafe bool CheckEncounter(Orbit p, Orbit nextPatch, double startEpoch, OrbitDriver sec, CelestialBody targetBody, PatchedConics.SolverParameters pars, bool logErrors)
        {
            try
            {
                Orbit s = sec.orbit;

                bool alwaysShowMarkers = GameSettings.ALWAYS_SHOW_TARGET_APPROACH_MARKERS && sec.celestialBody == targetBody;

                /*
                 * --- Stage 1: Pe/Ap prefilter ---
                 * Quick rejection based on whether the radial extent of the two orbits overlap, unless we are always showing target markers.
                 */

                double SoIbuffer = 1.1;

                if (!alwaysShowMarkers && !Orbit.PeApIntersects(p, s, sec.celestialBody.sphereOfInfluence * SoIbuffer))
                    return false;

                // these variables are almost entirely useless, and seem very useless now that this will trigger if alwaysShowMarkers is true
                if (p.closestEncounterLevel < Orbit.EncounterSolutionLevel.ORBIT_INTERSECT)
                {
                    p.closestEncounterLevel = Orbit.EncounterSolutionLevel.ORBIT_INTERSECT;
                    p.closestEncounterBody  = sec.celestialBody;
                }

                /*
                 * --- Stage 2: Geometric MOID solver ---
                 * This is now Baluev's solver and returns up to 4 minima (which should be as many as have ever been found?)
                 */

                // a = semi-major axis, e = eccentricity, i = inclination,
                // w = argument of periapsis, Om = longitude of ascending node
                var o1 = new COrbitData(p.semiMajorAxis, p.eccentricity, Deg2Rad(p.inclination), Deg2Rad(p.argumentOfPeriapsis), Deg2Rad(p.LAN));
                var o2 = new COrbitData(s.semiMajorAxis, s.eccentricity, Deg2Rad(s.inclination), Deg2Rad(s.argumentOfPeriapsis), Deg2Rad(s.LAN));

                Baluev.MoidInfo* info = stackalloc Baluev.MoidInfo[4];

                int num = Baluev.FindAllMinima(o1, o2, info);

                /*
                for (int i = 0; i < 4; i++)
                    Logger.Print($"{num} {info[i].dst} {info[i].u1} {info[i].u2}"); */

                // perfectly circular coplanar orbit check.
                if (num < 1)
                {
                    if (logErrors && !Thread.CurrentThread.IsBackground)
                        Debug.Log("CheckEncounter: failed to find any intercepts at all");

                    return false;
                }

                // resolve the true anomalies for the primary
                for (int i = 0; i < num; i++)
                    info[i].ta1 = p.GetTrueAnomaly(info[i].u1);

                // calculate the transit time
                for (int i = 0; i < num; i++)
                    info[i].tt = p.GetDTforTrueAnomaly(info[i].ta1, 0.0);

                // fix the transit time for possible bugs in GetDTforTrueAnomaly()
                if (p.eccentricity < 1.0)
                    for (int i = 0; i < num; i++)
                    {
                        double period = p.period;
                        info[i].tt -= period * Math.Floor(info[i].tt / period);
                    }

                // insertion sort by transit time
                for (int i = 1; i < num; i++)
                {
                    Baluev.MoidInfo keyInfo = info[i];

                    int j = i - 1;
                    while (j >= 0 && info[j].tt > keyInfo.tt)
                    {
                        info[j + 1] = info[j];
                        j--;
                    }

                    info[j + 1] = keyInfo;
                }

                // resolve true anomalies for secondary
                for (int i = 0; i < num; i++)
                    info[i].ta2 = s.GetTrueAnomaly(info[i].u2);

                // resolve UT of geometric point
                for (int i = 0; i < num; i++)
                    info[i].ut = startEpoch + info[i].tt;

                double bestClAppr = double.PositiveInfinity;
                double bestUTappr = 0;

                for (int i = 0; i < num; i++)
                {
                    if ((!alwaysShowMarkers && info[i].dst >= sec.celestialBody.sphereOfInfluence) || double.IsInfinity(info[i].ut))
                        continue;
                    if (CheckGeometricalEncounter(p, nextPatch, startEpoch, sec, pars, info[i], s, ref bestClAppr, ref bestUTappr))
                        return true;
                }

                if (p.eccentricity < 1.0)
                {
                    // try wrapping around the end of the orbit to try the first MOID with a seed one period in the future.
                    Baluev.MoidInfo wrappedInfo = info[0];
                    wrappedInfo.tt += p.period;
                    wrappedInfo.ut =  startEpoch + wrappedInfo.tt;

                    if ((alwaysShowMarkers || !(wrappedInfo.dst >= sec.celestialBody.sphereOfInfluence)) && !double.IsInfinity(wrappedInfo.ut))
                    {
                        if (CheckGeometricalEncounter(p, nextPatch, startEpoch, sec, pars, wrappedInfo, s, ref bestClAppr, ref bestUTappr))
                            return true;
                    }
                }

                if (bestClAppr < double.PositiveInfinity)
                    p.closestTgtApprUT = bestUTappr;

                return false;
            }

            catch (Exception value)
            {
                //Logger.Print($"{value}"); // TODO: remove this
                if (!Thread.CurrentThread.IsBackground)
                {
                    Console.WriteLine(value);
                }

                return false;
            }
        }

        private static bool CheckGeometricalEncounter(Orbit p, Orbit nextPatch, double startEpoch, OrbitDriver sec, PatchedConics.SolverParameters pars, Baluev.MoidInfo info, Orbit s, ref double bestClAppr, ref double bestUTappr)
        {
            // These APIs seem to be almost entirely unused throughout the stock codebase.  We never set
            // SOI_INTERSECT_2 here because that is tightly coupled to the idea of only have 2 possible minima.
            if (p.closestEncounterLevel < Orbit.EncounterSolutionLevel.SOI_INTERSECT_1)
            {
                p.closestEncounterLevel = Orbit.EncounterSolutionLevel.SOI_INTERSECT_1;
                p.closestEncounterBody  = sec.celestialBody;
            }

            // estimate how long it takes to cross the SOI
            double vrel  = (p.getOrbitalVelocityAtTrueAnomaly(info.ta1) - s.getOrbitalVelocityAtTrueAnomaly(info.ta2)).magnitude;
            double maxDT = sec.celestialBody.sphereOfInfluence / vrel;

            p.UTappr = info.ut;
            p.ClAppr = GetClosestApproach(p, s, startEpoch, maxDT, pars);

            // These APIs seem entire unused by stock code and tightly coupled to the idea of only having 2 minima points, so just
            // consistently fill them with the best approach data.  If this is ever found to be buggy with some mod we can figure out
            // what it expects these to have when there's 2, 3 or 4 minima.
            if (p.ClAppr < bestClAppr)
            {
                p.timeToTransition1         = info.tt;
                p.secondaryPosAtTransition1 = s.getPositionAtUT(info.ut);
                p.timeToTransition2         = info.tt;
                p.secondaryPosAtTransition2 = p.secondaryPosAtTransition1;
                p.nearestTT                 = info.tt;
                p.nextTT                    = info.tt;
                bestClAppr                  = p.ClAppr;
                bestUTappr                  = p.UTappr;
            }

            if (EncountersBody(p, s, nextPatch, sec, startEpoch, pars))
                return true;

            return false;
        }

        // This is the CheckEncounter solver running against the old MOID Solver, but with considerable bugfixing.  It is close
        // to the best accuracy the old MOID solver could get, but it still has bugs in the prefiltering based on the MOID time.
        public static bool CheckEncounterOldMOID(Orbit p, Orbit nextPatch, double startEpoch, OrbitDriver sec, CelestialBody targetBody, PatchedConics.SolverParameters pars, bool logErrors)
        {
            try
            {
                Orbit s = sec.orbit;

                bool alwaysShowMarkers = GameSettings.ALWAYS_SHOW_TARGET_APPROACH_MARKERS && sec.celestialBody == targetBody;

                double SoIbuffer = 1.1;

                if (!alwaysShowMarkers && !Orbit.PeApIntersects(p, s, sec.celestialBody.sphereOfInfluence * SoIbuffer))
                    return false;

                if (p.closestEncounterLevel < Orbit.EncounterSolutionLevel.ORBIT_INTERSECT)
                {
                    p.closestEncounterLevel = Orbit.EncounterSolutionLevel.ORBIT_INTERSECT;
                    p.closestEncounterBody  = sec.celestialBody;
                }

                double ClEctr1 = p.ClEctr1;
                double ClEctr2 = p.ClEctr2;
                double FEVp    = p.FEVp;
                double FEVs    = p.FEVs;
                double SEVp    = p.SEVp;
                double SEVs    = p.SEVs;

                int num = Orbit.FindClosestPoints(p, s, ref ClEctr1, ref ClEctr2, ref FEVp, ref FEVs, ref SEVp, ref SEVs, 0.0001, pars.maxGeometrySolverIterations, ref pars.GeoSolverIterations);

                //Logger.Print($"ClEctr1: {ClEctr1} ClEctr2: {ClEctr2} FEVp: {FEVp} FEVs: {FEVs} SEVp: {SEVp} SEVs: {SEVs}");

                if (num < 1)
                {
                    if (logErrors && !Thread.CurrentThread.IsBackground)
                        Debug.Log("CheckEncounter: failed to find any intercepts at all");

                    return false;
                }

                double tt1 = p.GetDTforTrueAnomaly(FEVp, 0.0);
                double tt2 = p.GetDTforTrueAnomaly(SEVp, 0.0);

                if (p.eccentricity < 1.0)
                {
                    double period = p.period;
                    tt1 -= period * Math.Floor(tt1 / period);
                    tt2 -= period * Math.Floor(tt2 / period);
                }

                double ut1 = tt1 + startEpoch;
                double ut2 = tt2 + startEpoch;

                if (double.IsInfinity(ut1) && double.IsInfinity(ut2))
                {
                    if (logErrors && !Thread.CurrentThread.IsBackground)
                        Debug.Log("CheckEncounter: both intercept UTs are infinite");

                    return false;
                }

                // XXX: This version still does known buggy pre-filtering based on the geometric time.
                if ((ut1 < p.StartUT || ut1 > p.EndUT) && (ut2 < p.StartUT || ut2 > p.EndUT))
                    return false;

                if (ut2 < ut1 || ut1 < p.StartUT || ut1 > p.EndUT)
                {
                    UtilMath.SwapValues(ref FEVp, ref SEVp);
                    UtilMath.SwapValues(ref FEVs, ref SEVs);
                    UtilMath.SwapValues(ref ClEctr1, ref ClEctr2);
                    UtilMath.SwapValues(ref tt1, ref tt2);
                    UtilMath.SwapValues(ref ut1, ref ut2);
                }

                if (ut2 < p.StartUT || ut2 > p.EndUT || double.IsInfinity(ut2))
                    num = 1;

                p.numClosePoints = num;
                p.FEVp           = FEVp;
                p.FEVs           = FEVs;
                p.SEVp           = SEVp;
                p.SEVs           = SEVs;
                p.ClEctr1        = ClEctr1;
                p.ClEctr2        = ClEctr2;

                if (p.closestEncounterLevel < Orbit.EncounterSolutionLevel.SOI_INTERSECT_1)
                {
                    p.closestEncounterLevel = Orbit.EncounterSolutionLevel.SOI_INTERSECT_1;
                    p.closestEncounterBody  = sec.celestialBody;
                }

                p.timeToTransition1         = tt1;
                p.secondaryPosAtTransition1 = s.getPositionAtUT(ut1);
                p.timeToTransition2         = tt2;
                p.secondaryPosAtTransition2 = s.getPositionAtUT(ut2);
                p.nearestTT                 = p.timeToTransition1;
                p.nextTT                    = p.timeToTransition2;

                if (double.IsNaN(p.nearestTT) && logErrors && !Thread.CurrentThread.IsBackground)
                    Debug.Log("nearestTT is NaN! t1: " + p.timeToTransition1 + ", t2: " + p.timeToTransition2 + ", FEVp: " + p.FEVp + ", SEVp: " + p.SEVp);

                double bestClAppr = double.PositiveInfinity;
                double bestUTappr = 0;

                if (alwaysShowMarkers || p.ClEctr1 < sec.celestialBody.sphereOfInfluence)
                {
                    double fVrel = (p.getOrbitalVelocityAtTrueAnomaly(FEVp)
                        - s.getOrbitalVelocityAtTrueAnomaly(FEVs)).magnitude;
                    double fMaxDT = sec.celestialBody.sphereOfInfluence / fVrel;

                    p.UTappr = startEpoch + p.nearestTT;
                    p.ClAppr = GetClosestApproach(p, s, startEpoch, 0.5 * fMaxDT, pars);

                    if (EncountersBody(p, s, nextPatch, sec, startEpoch, pars))
                        return true;

                    bestClAppr = p.ClAppr;
                    bestUTappr = p.UTappr;
                }

                if (num > 1)
                {
                    if (alwaysShowMarkers || p.ClEctr2 < sec.celestialBody.sphereOfInfluence)
                    {
                        p.closestEncounterLevel = Orbit.EncounterSolutionLevel.SOI_INTERSECT_2;
                        p.closestEncounterBody  = sec.celestialBody;

                        double sVrel = (p.getOrbitalVelocityAtTrueAnomaly(SEVp)
                            - s.getOrbitalVelocityAtTrueAnomaly(SEVs)).magnitude;
                        double sMaxDT = sec.celestialBody.sphereOfInfluence / sVrel;

                        p.UTappr = startEpoch + p.nextTT;
                        p.ClAppr = GetClosestApproach(p, s, startEpoch, 0.5 * sMaxDT, pars);

                        if (EncountersBody(p, s, nextPatch, sec, startEpoch, pars))
                            return true;

                        if (p.ClAppr < bestClAppr)
                        {
                            bestClAppr = p.ClAppr;
                            bestUTappr = p.UTappr;
                        }
                    }
                }

                if (bestClAppr < double.PositiveInfinity)
                    p.closestTgtApprUT = bestUTappr;

                return false;
            }

            catch (Exception value)
            {
                if (!Thread.CurrentThread.IsBackground)
                {
                    Console.WriteLine(value);
                }

                return false;
            }
        }

        /// <summary>
        ///     Wrapper function around SolveClosestApproach() to determine the MinUT/MaxUT bounds for the
        ///     the search.  This ultimately refines the geometric closest approach of the conic sections to
        ///     the actual time domain closest approach of the two orbits.
        ///     The search window is determined by the primary orbit's type:
        ///     - Elliptical: one full orbital period is searched [startEpoch, startEpoch + period].
        ///     - Hyperbolic: if the orbit is bounded by exiting an SOI boundary, use the time of that SOI
        ///     crossing as the upper limit.  If the parent is a root body without an SOI, then if the secondary
        ///     is elliptical (Planet) cap the search to 3xSMA of the Moon.  If both orbits are hyperbolic
        ///     ("Comet") then use an arbitrary upper limit based on the size of the solar system.
        ///     The result is stored in p.UTappr (the UT of closest approach) and Clappr value is returned (the
        ///     closest approach distance). These are used downstream by EncountersBody() to find the SOI crossing
        ///     given the closest time domain approach.
        /// </summary>
        /// <param name="p">The primary orbit (vessel orbit/patch being tested).</param>
        /// <param name="s">The secondary orbit (celestial body being tested for encounter).</param>
        /// <param name="startEpoch">The UT at the start of this patch — the earliest time to search.</param>
        /// <param name="maxDT">Maximum clamp on the timestep that the Halley solver can take.</param>
        /// <param name="pars">Solver parameters including iteration limits and convergence epsilon.</param>
        /// <returns>The closest approach distance found within the search window.</returns>
        public static double GetClosestApproach(Orbit p, Orbit s, double startEpoch, double maxDT, PatchedConics.SolverParameters pars)
        {
            double maxUT;

            if (p.eccentricity < 1.0)
            {
                // Elliptical orbit: search one full period of the orbit.
                maxUT = startEpoch + p.period;
            }
            else
            {
                // Vessel is on a hyperbolic orbit
                double taForSOI;

                // Handle the case where the reference body is the root (e.g. Kerbin/Sun)
                if (double.IsInfinity(p.referenceBody.sphereOfInfluence))
                {
                    if (s.eccentricity < 1.0)
                        // Transfer to a closed orbit (a Planet)
                        // XXX: this function doesn't have access to the planet's SOI, just its orbit, so use 3*SMA rather than Apoapsis+SOI
                        // (this should not be buggy and not important enough to change the method signature to fix)
                        taForSOI = p.TrueAnomalyAtRadius(s.semiMajorAxis * 3.0);
                    else
                        // Both orbits are hyperbolic (transfer to a "Comet"):  Just use a very large distance.
                        taForSOI = p.TrueAnomalyAtRadius(pars.outerReaches);
                }
                else
                {
                    // Ejection out of a Planet/Moon SOI: use the SOI transition as the natural end point of the search.
                    taForSOI = p.TrueAnomalyAtRadius(p.referenceBody.sphereOfInfluence);
                }

                // Note that by accident or design TrueAnomalyAtRadius() always gives positive TA, so picks the departure
                // hyperbola correctly, so once we're past this point, we're leaving forever...
                double utForSOI = p.GetUTforTrueAnomaly(taForSOI, 0.0);
                maxUT = utForSOI;
            }

            // Refine the geometric closest approach into the time domain closest approach using the bounds
            return SolveClosestApproach(p, s, ref p.UTappr, maxDT, 0.0, startEpoch, maxUT, pars.epsilon, pars.maxTimeSolverIterations, ref pars.TimeSolverIterations1);
        }

        static (double r, double rdv) GetRdvAtUT(Orbit p, Orbit s, double ut)
        {
            p.GetOrbitalStateVectorsAtUT(ut, out Vector3d pPos, out Vector3d pVel);
            s.GetOrbitalStateVectorsAtUT(ut, out Vector3d sPos, out Vector3d sVel);
            return ((pPos - sPos).magnitude, Vector3d.Dot(pVel - sVel, pPos - sPos));
        }

        /// <summary>
        ///     Finds the time of closest approach between two orbiting bodies using Halley's method
        ///     on the range-rate function.
        ///     The solver minimizes the squared separation distance |r_p(t) - r_s(t)|² by finding
        ///     the root of its time derivative, the range-rate function:
        ///     rdv   = d(|r|²/2)/dt   = dot(r_rel, v_rel)         — zero at closest approach
        ///     drdv  = d²(|r|²/2)/dt² = |v_rel|² + dot(r_rel, a_rel) — positive at a minimum
        ///     ddrdv = d³(|r|²/2)/dt³ = 3·dot(v_rel, a_rel) + dot(r_rel, j_rel)
        ///     The algorithm proceeds in three phases:
        ///     Phase 1 (Presolve): Advance UT into a basin of convergence around a minimum (where
        ///     drdv > 0). Uses the quadratic approximation of rdv to estimate zero-crossing steps,
        ///     with clamping to prevent overshooting. This phase handles the case where the initial
        ///     UT is far from any minimum — e.g. near a maximum of the distance function, or on
        ///     the wrong side of a minimum where the curvature is concave-down.
        ///     Phase 2 (Halley iteration): Once in the right basin, applies Halley's method for
        ///     cubic convergence. Each iteration roughly triples the number of correct digits.
        ///     Typically converges in 3–5 iterations.
        ///     Boundary handling: If the Halley step would exit [MinUT, MaxUT], the solver clamps
        ///     or wraps. For elliptical orbits where the search window spans a full period, a
        ///     single MaxDT bump is applied to allow the solver to wrap around the orbit.
        ///     IMPORTANT: This solver is a local minimizer — it converges to the nearest minimum from
        ///     the starting UT. It does not attempt to choose between multiple geometric intercept
        ///     points. That responsibility belongs to the caller (_CheckEncounter), which invokes
        ///     this solver once per candidate intercept with an appropriate initial UT.
        /// </summary>
        /// <param name="p">Primary orbit (vessel).</param>
        /// <param name="s">Secondary orbit (celestial body).</param>
        /// <param name="ut">
        ///     On entry, the initial guess for the time of closest approach.
        ///     On exit, the converged UT of the closest approach found.
        /// </param>
        /// <param name="maxDT">Maximum clamp on the timestep that the Halley solver can take.</param>
        /// <param name="threshold">Unused (from previous bisection implementation)</param>
        /// <param name="minUT">Earliest allowed UT (start of the current patch).</param>
        /// <param name="maxUT">Latest allowed UT (end of the search window).</param>
        /// <param name="epsilon">Convergence threshold — iteration stops when |dT| &lt; epsilon.</param>
        /// <param name="maxIterations">Hard iteration cap across all phases combined.</param>
        /// <param name="iterationCount">On exit, the total number of state evaluations performed.</param>
        /// <returns>
        ///     The closest approach distance (magnitude of relative position)
        /// </returns>
        public static double SolveClosestApproach(Orbit p, Orbit s, ref double ut, double maxDT, double threshold, double minUT, double maxUT, double epsilon, int maxIterations, ref int iterationCount)
        {
            /*
             * Establish a bracket
             */

            iterationCount = 0;

            (double _, double rdv) = GetRdvAtUT(p, s, ut);
            int startSign = Math.Sign(rdv);

            double step = maxDT / 8;

            iterationCount = 0;

            double bracketHi = double.NaN;
            double bracketLo = double.NaN;

            double probe = ut, probe1 = ut, probe2 = ut;
            while (true)
            {
                probe2   = probe1;
                probe1   = probe;
                probe    = ut + step;
                (_, rdv) = GetRdvAtUT(p, s, probe);
                if (rdv * startSign <= 0)
                {
                    if (rdv > 0)
                    {
                        bracketHi = probe;
                        bracketLo = probe2;
                    }
                    else
                    {
                        bracketLo = probe;
                        bracketHi = probe2;
                    }

                    break;
                }

                probe2   = probe1;
                probe1   = probe;
                probe    = ut - step;
                (_, rdv) = GetRdvAtUT(p, s, probe);
                if (rdv * startSign <= 0)
                {
                    if (rdv > 0)
                    {
                        bracketHi = probe;
                        bracketLo = probe2;
                    }
                    else
                    {
                        bracketLo = probe;
                        bracketHi = probe2;
                    }

                    break;
                }

                step *= 1.6;
                if (iterationCount++ > maxIterations)
                    throw new Exception("could not find bracket"); // TODO: this cannot throw
            }

            /*
            var state2 = new Orbit.CASolutionState(p, s, maxDT);

            for (int i = 0; i <= 100; i++)
            {
                double probeX = minUT + i / 100.0 * (maxUT - minUT);
                state2.Update(probeX, ref iterationCount, true);
                Logger.Print($"{probeX}: {state2.rstate.pos.magnitude} {state2.targetAhead} {state2.rdv} {state2.drdv}");
                double nextProbe = minUT + (i+1) / 100.0 * (maxUT - minUT);
                if (ut > probeX && ut < nextProbe)
                {
                    state2.Update(ut, ref iterationCount, true);
                    Logger.Print($"{ut}: {state2.rstate.pos.magnitude} {state2.targetAhead} {state2.rdv} {state2.drdv} [MOID]");
                }
                if (bracketLo > probeX && bracketLo < nextProbe)
                {
                    state2.Update(bracketLo, ref iterationCount, true);
                    Logger.Print($"{bracketLo}: {state2.rstate.pos.magnitude} {state2.targetAhead} {state2.rdv} {state2.drdv} [BRACKETLO]");
                }
                if (bracketHi > probeX && bracketHi < nextProbe)
                {
                    state2.Update(bracketHi, ref iterationCount, true);
                    Logger.Print($"{bracketHi}: {state2.rstate.pos.magnitude} {state2.targetAhead} {state2.rdv} {state2.drdv} [BRACKETHI]");
                }
            }
            */

            ut = (bracketLo + bracketHi) * 0.5;

            iterationCount = 0;

            /*
             * Main Halley iteration
             */

            var state = new Orbit.CASolutionState(p, s, maxDT);

            state.Update(ut, ref iterationCount);

            double halleyDt = state.Halley_dt();
            double newtonDt = -state.rdv / state.drdv;
            double dt       = Math.Sign(halleyDt) != Math.Sign(newtonDt) ? newtonDt : halleyDt;

            while (iterationCount < maxIterations && !(Math.Abs(dt) <= epsilon))
            {
                double candidateUT = ut + dt;

                // fallback bisection
                if (double.IsNaN(candidateUT) || candidateUT <= bracketLo || candidateUT >= bracketHi)
                    candidateUT = (bracketLo + bracketHi) * 0.5;

                ut = candidateUT;
                state.Update(ut, ref iterationCount, true);
                //Logger.Print($"HALLEY: {ut} halley:{halleyDt} newton:{newtonDt} targetAhead: {state.targetAhead} rdv:{state.rdv} drdv:{state.drdv}");

                if (state.rdv < 0.0 && double.IsNaN(bracketLo))
                    bracketLo = ut;
                else if (state.rdv > 0.0 && double.IsNaN(bracketHi))
                    bracketHi = ut;
                else
                    break; // landed exactly on the root

                halleyDt = state.Halley_dt();
                newtonDt = -state.rdv / state.drdv;
                dt       = Math.Sign(halleyDt) != Math.Sign(newtonDt) ? newtonDt : halleyDt;
            }

            if (iterationCount >= maxIterations && GameSettings.VERBOSE_DEBUG_LOG && !Thread.CurrentThread.IsBackground)
            {
                Debug.Log("[Orbit] SolveClosestApproach: solve took too many iterations, result incorrect");
            }

            if (iterationCount >= maxIterations)
                throw new Exception("too many iterations"); // TODO: cannot throw, delete this

            return state.rstate.pos.magnitude;
        }

        // EncountersBody — PatchedConics.cs
        //
        // Called after GetClosestApproach() has already refined p.UTappr and p.ClAppr for
        // the time domain.  This function decides whether that closest approach
        // actually constitutes an SOI entry, and if so it pins down the exact SOI-crossing
        // time, validates the approach direction, and commits the patch transition.
        //
        // Parameters:
        //   p          — the spacecraft orbit being propagated (the "patch" under evaluation)
        //   s          — the orbit of the candidate body (sec.orbit)
        //   nextPatch  — the orbit object to be initialised if an encounter is confirmed
        //   sec        — the OrbitDriver wrapping the candidate celestial body
        //   startEpoch — the UT at which this patch begins
        //   pars       — solver tuning parameters (iteration limits, epsilon, etc.)
        //
        // Returns true  → encounter confirmed; p and nextPatch have been fully committed.
        // Returns false → no encounter this window.

        public static bool EncountersBody(Orbit p, Orbit s, Orbit nextPatch, OrbitDriver sec, double startEpoch, PatchedConics.SolverParameters pars)
        {
            // A sentinel value of -1 means SolveClosestApproach failed.
            if (p.ClAppr < 0)
                return false;

            // -----------------------------------------------------------------------
            // GATE CHECK: is the closest approach distance inside the body's SoI?
            //
            // p.ClAppr is the scalar distance between the spacecraft and the body at
            // p.UTappr, as computed by the preceding GetClosestApproach call.
            // -----------------------------------------------------------------------
            if (p.ClAppr < sec.celestialBody.sphereOfInfluence)
            {
                // -------------------------------------------------------------------
                // PHASE 1 — Bisect to the exact SOI crossing time.
                //
                // p.UTappr is the time of closest approach, which is already *inside*
                // the SoI.  We want the earlier moment when the spacecraft first
                // crossed the SoI boundary.
                //
                // SolveSOI uses binary-search / bisection (BSP variant removed in
                // -newer; plain SolveSOI used instead) over the interval
                // [startEpoch, p.UTappr].  The search half-width seed is
                // (p.UTappr - startEpoch) * 0.5, i.e. the midpoint of that interval,
                // which is a reasonable starting bracket.
                //
                // Result is written into p.UTsoi.
                // -------------------------------------------------------------------
                p.UTsoi = p.UTappr; // initialise bisection anchor at closest approach
                SolveSOI(
                    p, s,
                    ref p.UTsoi,                         // in/out: SOI crossing time
                    (p.UTappr - startEpoch) * 0.5,       // initial search half-width
                    sec.celestialBody.sphereOfInfluence, // target radius
                    startEpoch,                          // search lower bound
                    p.UTappr,                            // search upper bound
                    pars.epsilon,                        // convergence tolerance
                    pars.maxTimeSolverIterations,
                    ref pars.TimeSolverIterations1); // iteration counter (diagnostics)

                // -------------------------------------------------------------------
                // PHASE 2 — Validate approach direction at the SOI crossing time.
                //
                // SolveSOI can in principle converge on an *outbound* SOI crossing
                // (spacecraft leaving the body's vicinity) rather than an inbound one.
                // We reject that case by checking the sign of the radial velocity in
                // the body's reference frame.
                //
                // Relative position vector:  lhs = pos_spacecraft − pos_body
                // Relative velocity vector:  rhs = vel_spacecraft − vel_body
                //
                // Dot(lhs, rhs) > 0  →  spacecraft is moving away from the body
                //                        (outbound crossing — reject)
                // Dot(lhs, rhs) < 0  →  spacecraft is approaching the body
                //                        (inbound crossing — accept)
                // Dot(lhs, rhs) = 0  →  tangential; conservative choice is to reject
                // -------------------------------------------------------------------
                p.GetOrbitalStateVectorsAtUT(p.UTsoi, out Vector3d pos, out Vector3d vel);
                s.GetOrbitalStateVectorsAtUT(p.UTsoi, out Vector3d pos2, out Vector3d vel2);

                Vector3d lhs = pos - pos2; // relative position (spacecraft w.r.t. body)
                Vector3d rhs = vel - vel2; // relative velocity

                if (Vector3d.Dot(lhs, rhs) >= 0.0)
                {
                    // Moving apart at the alleged SOI crossing — this is an outbound
                    // event (or we started inside the SoI on an escaping trajectory).
                    // Discard; the caller may try the second intercept window.
                    return false;
                }

                // -------------------------------------------------------------------
                // PHASE 3 — Commit the encounter.
                //
                // Everything checks out: the spacecraft enters the body's SoI at
                // p.UTsoi on an inbound trajectory.  Initialise the next patch to
                // represent the orbit in the body's reference frame starting at that
                // moment, then seal off the current patch at the SOI boundary.
                // -------------------------------------------------------------------

                // Build the next patch: orbit in sec's reference frame at p.UTsoi.
                nextPatch.UpdateFromOrbitAtUT(p, p.UTsoi, sec.celestialBody);

                // Seal the current patch: it runs from its start epoch to the SOI entry.
                p.StartUT = startEpoch;
                p.EndUT   = p.UTsoi;

                // Mark the transition type so the solver knows why this patch ended.
                p.patchEndTransition = Orbit.PatchTransitionType.ENCOUNTER;

                return true;
            }

            // -----------------------------------------------------------------------
            // Closest approach is outside the SoI (or ClAppr == -1 sentinel).
            // No encounter this window.
            // -----------------------------------------------------------------------
            return false;
        }

        public static bool SolveSOI(Orbit p, Orbit s, ref double ut, double dT, double rsoi, double minUT, double maxUT, double epsilon, int maxIterations, ref int iterationCount)
        {
            double rdv       = Orbit.RelativeStateAtUT(p, s, maxUT, out Orbit.State _, out Orbit.State _, out Orbit.State rstate);
            double magnitude = rstate.pos.magnitude;

            iterationCount = 1;
            // ReSharper disable once CompareOfFloatsByEqualityOperator
            if (magnitude == rsoi)
            {
                ut = maxUT;
                return true;
            }

            UtilMath.SphereIntersection(rsoi, rstate.pos, rstate.vel, out dT, false);
            ut = Math.Max(minUT, maxUT + dT);

            double rsoi2  = rsoi * rsoi;

            while (iterationCount++ < maxIterations)
            {
                rdv    = Orbit.RelativeStateAtUT(p, s, ut, out _, out _, out rstate);
                double sqrMag = rstate.pos.sqrMagnitude;

                if (NearlyEqual(sqrMag, rsoi2, 1e-10))
                    break;

                if (sqrMag < rsoi2)
                    maxUT = ut;
                else
                    minUT = ut;

                dT = (rsoi2 - sqrMag) / (2.0 * rdv);
                double next = ut + dT;

                if (next <= minUT || next >= maxUT)
                    next = 0.5 * (minUT + maxUT);

                // ReSharper disable once CompareOfFloatsByEqualityOperator
                if (next == ut)
                    break;

                ut = next;
            }

            return iterationCount < maxIterations;
        }
    }
}
