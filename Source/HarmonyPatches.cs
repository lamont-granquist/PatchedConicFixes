using System;
using System.Threading;
using HarmonyLib;
using UnityEngine;

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
            __result = HarmonyPatches.CheckEncounter(p, nextPatch, startEpoch, sec, targetBody, pars, logErrors);

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

    public static class HarmonyPatches
    {
        public static string OrbitString(Orbit o) => $"{o.inclination}, {o.eccentricity}, {o.semiMajorAxis}, {o.LAN}, {o.argumentOfPeriapsis}, {o.meanAnomalyAtEpoch}, {o.epoch}";

        /// <summary>
        ///     Finds the closest approach distance between two orbits within a bounded time window, using a
        ///     BSP (Binary Space Partitioning) search in the time domain.
        ///     The search window is determined by the primary orbit's type:
        ///     - Elliptical (e &lt; 1): The entire orbital period is searched, since the vessel will revisit
        ///     every point on the orbit. The window spans [startEpoch, startEpoch + period].
        ///     - Hyperbolic (e ≥ 1): The orbit is open and extends to infinity, but the BSP solver cannot
        ///     subdivide an infinite interval (inf/2 = inf, inf - inf = NaN). Instead, the search is
        ///     bounded by computing the true anomaly at which the vessel reaches the reference body's SOI
        ///     boundary, and converting that to a UT. Beyond the SOI the vessel will be on a different
        ///     patch anyway, so searching further is pointless.
        ///     - For hyperbolic orbits around the root body (infinite SOI), we use heuristic bounds: 3x the
        ///     secondary body's SMA if it has a closed orbit, or pars.outerReaches as a last resort.
        ///     The result is stored in p.UTappr (the UT of closest approach) and returned as p.ClAppr (the
        ///     closest approach distance). These are used downstream by _EncountersBody to determine whether
        ///     the vessel enters the secondary body's SOI.
        /// </summary>
        /// <param name="p">The primary orbit (vessel or patch being tested).</param>
        /// <param name="s">The secondary orbit (celestial body being tested for encounter).</param>
        /// <param name="startEpoch">The UT at the start of this patch — the earliest time to search.</param>
        /// <param name="dT">
        ///     Initial step size for the BSP solver. Passed through from the caller, typically
        ///     half the time-to-transition to the nearest geometric intercept point.
        /// </param>
        /// <param name="pars">Solver parameters including iteration limits and convergence epsilon.</param>
        /// <returns>The closest approach distance found within the search window.</returns>
        public static double GetClosestApproach(Orbit p, Orbit s, double startEpoch, double dT, PatchedConics.SolverParameters pars)
        {
            double MinUT;
            double MaxUT;

            if (p.eccentricity < 1.0)
            {
                // Elliptical orbit: the vessel will traverse the entire orbit within one period,
                // so search the full period starting from the patch epoch. This guarantees we don't
                // miss an encounter that occurs late in the orbit.
                MinUT = startEpoch;
                MaxUT = startEpoch + p.period;
            }
            else
            {
                // Hyperbolic/parabolic orbit: the trajectory extends to infinity, but the BSP solver
                // requires finite bounds. Compute the true anomaly at which the vessel reaches the
                // outermost physically meaningful distance, then convert to UT for the search bound.
                double AforSOI;

                if (double.IsInfinity(p.referenceBody.sphereOfInfluence))
                {
                    // The reference body is the root body (e.g. the Sun in stock KSP) which has an
                    // infinite SOI. We need a heuristic upper bound for the search.
                    if (s.eccentricity < 1.0)
                    {
                        // The secondary body has a closed orbit, so there's no point searching beyond
                        // roughly 3x its SMA — the vessel will be well past any possible encounter.
                        AforSOI = p.TrueAnomalyAtRadius(s.semiMajorAxis * 3.0);
                    }
                    else
                    {
                        // Both orbits are open (e.g. an escape trajectory past an escaping body).
                        // This case doesn't arise in stock KSP. Fall back to a very large distance.
                        AforSOI = p.TrueAnomalyAtRadius(pars.outerReaches);
                    }
                }
                else
                {
                    // Standard case: bound the search at the reference body's SOI. Beyond this the
                    // vessel transitions to a different patch, so the current orbit is no longer valid.
                    AforSOI = p.TrueAnomalyAtRadius(p.referenceBody.sphereOfInfluence);
                }

                double TforSOI = p.GetUTforTrueAnomaly(AforSOI, 0.0);
                MinUT = startEpoch;
                MaxUT = TforSOI;
            }

            // Run the BSP solver over the computed time window. SolveClosestApproach evaluates the
            // 3D separation |p.position(UT) - s.position(UT)| at candidate times and converges on
            // the minimum via ternary search / bisection. The result is written to p.UTappr (time)
            // and returned as the distance.
            return SolveClosestApproach(p, s, ref p.UTappr, dT, 0.0, MinUT, MaxUT, pars.epsilon, pars.maxTimeSolverIterations, ref pars.TimeSolverIterations1);
        }

        /// <summary>
        ///     Determines whether the primary orbit encounters (enters the SOI of) the given celestial body.
        ///     This is the main encounter detection entry point, implementing a multi-stage filter pipeline:
        ///     Stage 1 — Pe/Ap prefilter: Checks whether the altitude bands of the two orbits overlap
        ///     (with a buffer for the body's SOI radius). This is a cheap scalar comparison that eliminates
        ///     obviously non-intersecting orbit pairs. For the targeted body, the buffer is widened using
        ///     sqrt(SMA/SOI) so that approach markers can be shown even for near-miss trajectories.
        ///     Stage 2 — Geometric MOID solver: Calls FindClosestPoints to find the one or two points on
        ///     the orbits where the inter-orbit distance function has critical points (local minima). These
        ///     are purely geometric — they identify WHERE on the orbits a close approach could occur,
        ///     independent of when the bodies actually reach those points.
        ///     Stage 3 — Time validation: Converts the geometric true anomalies to universal times and
        ///     checks whether they fall within the current patch's time bounds [StartUT, EndUT]. Intercepts
        ///     outside the patch are discarded. The remaining intercepts are sorted chronologically.
        ///     Stage 4 — Time-domain closest approach: For each valid geometric intercept (up to two),
        ///     runs the time-domain BSP solver (GetClosestApproach) seeded near that intercept to find the
        ///     actual minimum separation accounting for both bodies' motion. If the closest approach is
        ///     within the body's SOI, calls EncountersBody to compute the SOI crossing time and generate
        ///     the next orbit patch.
        ///     IMPORTANT: Both intercept points must be checked. The geometric solver finds critical points
        ///     of the distance function between the two orbital paths, which typically yields two candidates
        ///     on roughly opposite sides of the orbits. Which candidate produces an actual SOI encounter
        ///     depends on the phasing of the bodies at the current epoch — information the geometry solver
        ///     doesn't have. If only the first intercept is checked, encounters near the second intercept
        ///     are silently missed, causing the well-known "phantom encounter flickering" bug where an
        ///     encounter appears and disappears as the player adjusts a maneuver node.
        /// </summary>
        /// <param name="p">
        ///     The primary orbit (vessel patch being tested). Modified in place with encounter
        ///     data if an encounter is found (UTappr, ClAppr, UTsoi, EndUT, patchEndTransition, etc).
        /// </param>
        /// <param name="nextPatch">
        ///     The orbit object to populate with the post-encounter trajectory if an
        ///     SOI transition occurs.
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
        public static bool CheckEncounter(Orbit p, Orbit nextPatch, double startEpoch, OrbitDriver sec, CelestialBody targetBody, PatchedConics.SolverParameters pars, bool logErrors)
        {
            try
            {
                Orbit s = sec.orbit;

                bool alwaysShowMarkers = GameSettings.ALWAYS_SHOW_TARGET_APPROACH_MARKERS && sec.celestialBody == targetBody;

                // --- Stage 1: Pe/Ap prefilter ---
                // Quick rejection based on whether the radial extent of the two orbits overlap.
                // If the primary's periapsis is above the secondary's apoapsis (plus SOI buffer),
                // or the primary's apoapsis is below the secondary's periapsis, there's no overlap
                // and no encounter is possible.
                //
                // For the targeted body, the buffer is widened using sqrt(SMA/SOI) to ensure approach
                // markers are displayed even for trajectories that come close without entering the SOI.
                // This scaling reflects that for bodies with large orbits relative to their SOI (e.g.
                // Jool), a wider search radius is needed to catch displayable near-misses.
                double SoIbuffer = 1.1;

                if (!alwaysShowMarkers && !Orbit.PeApIntersects(p, s, sec.celestialBody.sphereOfInfluence * SoIbuffer))
                    return false;

                // If we get past the Pe/Ap check, the orbits at least overlap radially — record this
                // as the current best encounter solution level if nothing better has been found yet.
                if (p.closestEncounterLevel < Orbit.EncounterSolutionLevel.ORBIT_INTERSECT)
                {
                    p.closestEncounterLevel = Orbit.EncounterSolutionLevel.ORBIT_INTERSECT;
                    p.closestEncounterBody  = sec.celestialBody;
                }

                // --- Stage 2: Geometric MOID solver ---
                // Find the critical points of the distance function between the two orbital paths.
                // This is a purely geometric calculation — it finds the true anomaly values on each
                // orbit where the inter-orbit distance is locally minimized (or maximized).
                //
                // The results are copied into local variables first to avoid overwriting a previously
                // found valid encounter's data. If this candidate turns out to be invalid, the patch's
                // existing encounter info is preserved.
                //
                // FindClosestPoints returns 1 or 2, indicating how many distinct close-approach
                // regions were found. Typically there are two, on roughly opposite sides of the orbits.
                // Returns 0 only for degenerate cases (coplanar circular orbits).
                double ClEctr1 = p.ClEctr1; // closest distance at first intercept point
                double ClEctr2 = p.ClEctr2; // closest distance at second intercept point
                double FEVp    = p.FEVp;    // true anomaly on primary orbit at first intercept
                double FEVs    = p.FEVs;    // true anomaly on secondary orbit at first intercept
                double SEVp    = p.SEVp;    // true anomaly on primary orbit at second intercept
                double SEVs    = p.SEVs;    // true anomaly on secondary orbit at second intercept

                int num = Orbit.FindClosestPoints(p, s, ref ClEctr1, ref ClEctr2, ref FEVp, ref FEVs, ref SEVp, ref SEVs, 0.0001, pars.maxGeometrySolverIterations, ref pars.GeoSolverIterations);

                if (num < 1)
                {
                    // No intercepts found — this should only happen for perfectly circular coplanar
                    // orbits where the distance function has no local minimum (it's constant).
                    if (logErrors && !Thread.CurrentThread.IsBackground)
                        Debug.Log("CheckEncounter: failed to find any intercepts at all");

                    return false;
                }

                /*
                Vector3d r0 = p.getRelativePositionFromTrueAnomaly(FEVp);
                Vector3d r1 = s.getRelativePositionFromTrueAnomaly(FEVs);
                Vector3d v0 = p.getOrbitalVelocityAtTrueAnomaly(FEVp);
                Vector3d v1 = s.getOrbitalVelocityAtTrueAnomaly(FEVs);

                Vector3d dist = r0 - r1;

                double dot0 = Vector3d.Dot(dist.normalized, v0.normalized);
                double dot1 = Vector3d.Dot(dist.normalized, v1.normalized);

                Logger.Print($"MOID check: reported={ClEctr1:F1} actual={dist.magnitude:F1} " +
                    $"dot0={dot0:F6} dot1={dot1:F6}");

                if (num > 1)
                {
                    Vector3d r02 = p.getRelativePositionFromTrueAnomaly(SEVp);
                    Vector3d r12 = s.getRelativePositionFromTrueAnomaly(SEVs);
                    Vector3d v02 = p.getOrbitalVelocityAtTrueAnomaly(SEVp);
                    Vector3d v12 = s.getOrbitalVelocityAtTrueAnomaly(SEVs);

                    Vector3d dist2 = r02 - r12;

                    double dot02 = Vector3d.Dot(dist2.normalized, v02.normalized);
                    double dot12 = Vector3d.Dot(dist2.normalized, v12.normalized);

                    Logger.Print($"MOID check2: reported={ClEctr2:F1} actual={dist2.magnitude:F1} " +
                        $"dot02={dot02:F6} dot12={dot12:F6}");
                }
                */

                // --- Stage 3: Time validation ---
                // Convert the geometric true anomalies to universal times so we can check whether the
                // intercepts actually fall within this patch's time bounds. A geometric intercept that
                // occurs outside [StartUT, EndUT] is unreachable on this patch and must be discarded.
                double tt1 = p.GetDTforTrueAnomaly(FEVp, 0.0); // delta-time from epoch to first intercept
                double tt2 = p.GetDTforTrueAnomaly(SEVp, 0.0); // delta-time from epoch to second intercept
                double ut1 = tt1 + startEpoch;                 // absolute UT of first intercept
                double ut2 = tt2 + startEpoch;                 // absolute UT of second intercept

                if (double.IsInfinity(ut1) && double.IsInfinity(ut2))
                {
                    // Both intercepts are at infinite time — can happen for degenerate geometries
                    // on hyperbolic orbits where the intercept true anomaly is beyond the asymptote.
                    if (logErrors && !Thread.CurrentThread.IsBackground)
                        Debug.Log("CheckEncounter: both intercept UTs are infinite");

                    return false;
                }

                // Reject if neither intercept falls within the patch's time bounds.
                if ((ut1 < p.StartUT || ut1 > p.EndUT) && (ut2 < p.StartUT || ut2 > p.EndUT))
                    return false;

                // Sort the two intercepts so that the first (FEV) is the earlier valid one. If the
                // "first" intercept is actually later or out of bounds, swap the two so downstream
                // code can process them in chronological order. This ensures we always test the
                // earlier encounter opportunity first.
                if (ut2 < ut1 || ut1 < p.StartUT || ut1 > p.EndUT)
                {
                    UtilMath.SwapValues(ref FEVp, ref SEVp);
                    UtilMath.SwapValues(ref FEVs, ref SEVs);
                    UtilMath.SwapValues(ref ClEctr1, ref ClEctr2);
                    UtilMath.SwapValues(ref tt1, ref tt2);
                    UtilMath.SwapValues(ref ut1, ref ut2);
                }

                // If the second intercept is out of bounds or infinite, reduce to single-intercept mode.
                if (ut2 < p.StartUT || ut2 > p.EndUT || double.IsInfinity(ut2))
                    num = 1;

                // At least one intercept is valid — commit the geometry results to the patch.
                p.numClosePoints = num;
                p.FEVp           = FEVp;
                p.FEVs           = FEVs;
                p.SEVp           = SEVp;
                p.SEVs           = SEVs;
                p.ClEctr1        = ClEctr1;
                p.ClEctr2        = ClEctr2;

                // Check if the geometric closest distances are already beyond the body's SOI. If so,
                // no encounter is possible, but we may still want to record the closest approach time
                // for display purposes (target approach markers).
                /*
                if (Math.Min(p.ClEctr1, p.ClEctr2) > sec.celestialBody.sphereOfInfluence)
                {
                    if (GameSettings.ALWAYS_SHOW_TARGET_APPROACH_MARKERS && sec.celestialBody == targetBody)
                    {
                        p.UTappr           = startEpoch;
                        p.ClAppr           = GetClosestApproach(p, s, startEpoch, p.nearestTT * 0.5, pars);
                        p.closestTgtApprUT = p.UTappr;
                    }

                    return false;
                }
                */

                // The geometric distances suggest an SOI encounter is plausible — upgrade the
                // encounter solution level.
                if (p.closestEncounterLevel < Orbit.EncounterSolutionLevel.SOI_INTERSECT_1)
                {
                    p.closestEncounterLevel = Orbit.EncounterSolutionLevel.SOI_INTERSECT_1;
                    p.closestEncounterBody  = sec.celestialBody;
                }

                // Record transition data for both intercepts (used by orbit rendering and UI).
                p.timeToTransition1         = tt1;
                p.secondaryPosAtTransition1 = s.getPositionAtUT(ut1);
                p.timeToTransition2         = tt2;
                p.secondaryPosAtTransition2 = s.getPositionAtUT(ut2);
                p.nearestTT                 = p.timeToTransition1;
                p.nextTT                    = p.timeToTransition2;

                if (double.IsNaN(p.nearestTT) && logErrors && !Thread.CurrentThread.IsBackground)
                    Debug.Log("nearestTT is NaN! t1: " + p.timeToTransition1 + ", t2: " + p.timeToTransition2 + ", FEVp: " + p.FEVp + ", SEVp: " + p.SEVp);

                // --- Stage 4: Time-domain closest approach and SOI check ---
                // For each valid geometric intercept, run the time-domain BSP solver to find the actual
                // minimum 3D separation (accounting for both bodies' orbital motion), then check whether
                // that closest approach falls within the body's SOI.
                //
                // The initial guess for the solver (p.UTappr) is set to startEpoch, and the step size
                // (dT) is set to half the time-to-transition. This seeds the solver in a region where
                // the first intercept's encounter should be found.

                // --- Try the first (earlier) intercept ---

                double bestClAppr = double.PositiveInfinity;
                double bestUTappr = 0;

                if (alwaysShowMarkers || p.ClEctr1 < sec.celestialBody.sphereOfInfluence)
                {
                    // IMPORTANT FIX:  compute maximum allowed step size based on velocities of the orbits at the
                    // geometric MOID intersection point, and the size of the SOI.
                    double fVrel = (p.getOrbitalVelocityAtTrueAnomaly(FEVp)
                        - s.getOrbitalVelocityAtTrueAnomaly(FEVs)).magnitude;
                    double fMaxDT = sec.celestialBody.sphereOfInfluence / fVrel;

                    // IMPORTANT FIX:  seed start time at exactly the geometric MOID point
                    p.UTappr = startEpoch + p.nearestTT;
                    p.ClAppr = GetClosestApproach(p, s, startEpoch, fMaxDT, pars);

                    if (EncountersBody(p, s, nextPatch, sec, startEpoch, pars))
                        return true;

                    // if there's no encounter record closest encounter
                    bestClAppr         = p.ClAppr;
                    bestUTappr         = p.UTappr;
                }

                // --- Try the second (later) intercept, if one exists ---
                // This is critical: the two geometric intercepts represent fundamentally different
                // regions of the orbit pair where an encounter can occur. The first intercept may fail
                // to produce an SOI encounter because the bodies are out of phase at that point in their
                // orbits, while the second intercept may succeed (or vice versa). Both must be checked.
                //
                // The solver is re-seeded with UTappr at the first intercept time, and the step size
                // spans the gap between the two intercepts, directing the search toward the second
                // encounter region.
                if (num > 1)
                {
                    if (alwaysShowMarkers || p.ClEctr2 < sec.celestialBody.sphereOfInfluence)
                    {
                        p.closestEncounterLevel = Orbit.EncounterSolutionLevel.SOI_INTERSECT_2;
                        p.closestEncounterBody  = sec.celestialBody;

                        // IMPORTANT FIX:  compute maximum allowed step size based on velocities of the orbits at the
                        // geometric MOID intersection point, and the size of the SOI.
                        double sVrel = (p.getOrbitalVelocityAtTrueAnomaly(SEVp)
                            - s.getOrbitalVelocityAtTrueAnomaly(SEVs)).magnitude;
                        double sMaxDT = sec.celestialBody.sphereOfInfluence / sVrel;

                        // IMPORTANT FIX:  seed start time at exactly the geometric MOID point
                        p.UTappr = startEpoch + p.nextTT;
                        p.ClAppr = GetClosestApproach(p, s, startEpoch, sMaxDT, pars);

                        if (EncountersBody(p, s, nextPatch, sec, startEpoch, pars))
                            return true;

                        // if there's no encounter record closest encounter
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
        /// <param name="UT">
        ///     On entry, the initial guess for the time of closest approach.
        ///     On exit, the converged UT of the closest approach found.
        /// </param>
        /// <param name="dT">
        ///     Initial step size for the presolve phase. Typically half the
        ///     time-to-transition from the geometric intercept.
        /// </param>
        /// <param name="threshold">Unused in current implementation (reserved).</param>
        /// <param name="MinUT">Earliest allowed UT (start of the current patch).</param>
        /// <param name="MaxUT">Latest allowed UT (end of the search window).</param>
        /// <param name="epsilon">Convergence threshold — iteration stops when |dT| &lt; epsilon.</param>
        /// <param name="maxIterations">Hard iteration cap across all phases combined.</param>
        /// <param name="iterationCount">On exit, the total number of state evaluations performed.</param>
        /// <returns>
        ///     The closest approach distance (magnitude of relative position), or -1.0 if the
        ///     initial UT is outside [MinUT, MaxUT].
        /// </returns>
        public static double SolveClosestApproach(Orbit p, Orbit s, ref double UT, double dT, double threshold, double MinUT, double MaxUT, double epsilon, int maxIterations, ref int iterationCount)
        {
            if (UT < MinUT)
            {
                return -1.0;
            }

            if (UT > MaxUT)
            {
                return -1.0;
            }

            var state = new Orbit.CASolutionState(p, s, dT);
            iterationCount = 0;

            // Evaluate the range-rate function and its derivatives at the initial guess.
            state.Update(UT, ref iterationCount);

            // Handle the case where the target body is ahead of us along the orbit (dot(r, v_ship) < 0)
            // but we're receding (rdv > 0) with increasing recession rate (drdv > 0). This means we're
            // on the far side of a distance maximum — the closest approach is behind us in time.
            // Take a Newton step backward on rdv to get into the neighborhood of the minimum.
            if (state.targetAhead && state.rdv > 0.0 && state.drdv > 0.0)
            {
                dT =  -state.rdv / state.drdv;
                dT =  Math.Max(dT, -state.maxdt);
                UT += dT;
                state.Update(UT, ref iterationCount);
            }

            // --- Phase 1: Presolve ---
            // Iterate until drdv > 0, which indicates we're in a basin of convergence around a
            // local minimum of the distance function (not a maximum or saddle point).
            //
            // At each step, we solve the quadratic approximation of rdv for its zero crossing:
            //   rdv(dt) ≈ rdv + drdv·dt + (ddrdv/2)·dt² = 0
            //
            // Rearranging: (ddrdv)·dt² + 2·drdv·dt + 2·rdv = 0
            //
            // The discriminant is: drdv² - 2·ddrdv·rdv
            // The two roots are: (-drdv ± sqrt(discriminant)) / ddrdv
            //
            // Which root to use depends on the sign of rdv:
            //   rdv > 0 (receding): we need to step toward the zero crossing where rdv transitions
            //     from positive to negative. The larger-magnitude root reaches the minimum.
            //   rdv ≤ 0 (approaching): we need the nearer zero crossing where rdv transitions from
            //     negative back to zero at the minimum. The smaller-magnitude root suffices.
            //
            // When the discriminant is negative (no real roots), the quadratic approximation shows
            // rdv doesn't cross zero nearby. Fall back to a clamped step in the direction that
            // should reduce |rdv|.
            while (!(state.drdv > 0.0) && iterationCount < maxIterations)
            {
                double rdv   = state.rdv;
                double drdv  = state.drdv;
                double ddrdv = state.ddrdv;

                double discriminant = drdv * drdv - 2.0 * ddrdv * rdv;

                if (discriminant >= 0.0)
                {
                    double sqrtDisc = Math.Sqrt(discriminant);

                    if (rdv > 0.0)
                    {
                        // Receding from the body: take the larger root, which steps us past
                        // the distance maximum and toward the minimum on the far side.
                        dT = (-drdv + sqrtDisc) / ddrdv;
                    }
                    else
                    {
                        // Approaching the body (or at minimum): take the smaller root, which
                        // steps us to the nearer zero crossing — the actual closest approach.
                        dT = (-drdv - sqrtDisc) / ddrdv;
                    }

                    dT = state.Clamp_dt(dT);
                }
                else
                {
                    // No real roots — rdv doesn't cross zero in the quadratic approximation.
                    // Step in the direction that should move us toward a region where it does:
                    // if approaching (rdv ≤ 0), step forward; if receding (rdv > 0), step backward.
                    dT = rdv > 0.0 ? -state.maxdt : state.maxdt;
                }

                UT += dT;
                state.Update(UT, ref iterationCount);
            }

            if (iterationCount >= maxIterations)
            {
                if (GameSettings.VERBOSE_DEBUG_LOG && !Thread.CurrentThread.IsBackground)
                {
                    Debug.LogFormat("[Orbit] SolveClosestApproach: presolve took too many iterations, bailing UT:{0} MinUT:{1} MaxUT:{2}", UT, MinUT, MaxUT);
                }

                return state.rstate.pos.magnitude;
            }

            // --- Phase 2: Halley iteration ---
            // We now have drdv > 0, so we're in the basin of a distance minimum. Apply Halley's
            // method to find the root of rdv with cubic convergence:
            //
            //   dt = -2·rdv·drdv / (2·drdv² - rdv·ddrdv)
            //
            // This is a third-order method that uses rdv, drdv, and ddrdv to converge much faster
            // than Newton's method. Typically reaches machine precision in 3–5 iterations.
            dT = state.Halley_dt();

            // Boundary handling: if the Halley step would go below MinUT, check whether we can
            // wrap around. For elliptical orbits the search window spans a full period, so adding
            // MaxDT effectively wraps to the other side of the orbit. For hyperbolic orbits where
            // there's no periodicity, clamp to MinUT and return.
            if (UT + dT < MinUT)
            {
                if (!(MaxUT - MinUT >= state.MaxDT) || (!(p.eccentricity < 1.0) && !(s.eccentricity < 1.0)))
                {
                    UT = MinUT;
                    state.Update(UT, ref iterationCount, true);
                    return state.rstate.pos.magnitude;
                }

                dT += state.MaxDT;
            }
            else if (UT + dT > MaxUT)
            {
                UT = MaxUT;
                state.Update(UT, ref iterationCount, true);
                return state.rstate.pos.magnitude;
            }

            // Main Halley loop: iterate until the step size falls below epsilon.
            while (iterationCount < maxIterations && !(Math.Abs(dT) <= epsilon))
            {
                UT += dT;
                UT =  Math.Min(MaxUT, Math.Max(MinUT, UT));
                state.Update(UT, ref iterationCount, true);
                dT = state.Halley_dt();
            }

            if (iterationCount >= maxIterations && GameSettings.VERBOSE_DEBUG_LOG && !Thread.CurrentThread.IsBackground)
            {
                Debug.Log("[Orbit] SolveClosestApproach: solve took too many iterations, result incorrect");
            }

            return state.rstate.pos.magnitude;
        }

        // _EncountersBody — PatchedConics.cs
        //
        // Called after GetClosestApproach has already refined p.UTappr and p.ClAppr for
        // a specific intercept window.  This function decides whether that closest approach
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

        internal static bool EncountersBody(Orbit p, Orbit s, Orbit nextPatch, OrbitDriver sec, double startEpoch, PatchedConics.SolverParameters pars)
        {
            // -----------------------------------------------------------------------
            // GATE CHECK: is the closest approach distance inside the body's SoI?
            //
            // p.ClAppr is the scalar distance between the spacecraft and the body at
            // p.UTappr, as computed by the preceding GetClosestApproach call.
            // The sentinel value -1 means SolveClosestApproach failed to converge or
            // was not run; we must not treat that as "zero distance".
            // -----------------------------------------------------------------------
            if (p.ClAppr < sec.celestialBody.sphereOfInfluence && p.ClAppr != -1.0)
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
                Orbit.SolveSOI(
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
    }
}
