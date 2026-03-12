# Changelog

## 2026-03-11 - EncountersBody() Bug Fix

Revert from `SolveSOI()` to `SolveSOI_BSP()`.  The Newton's method solver in `SolveSOI()` has no protection against the derivative being near zero in a step, and it fails to converge, which causes known reported issues with encounters "flickering" and this has been confirmed to really solve an existing bug.  Maybe the Newton solver could be rescued, but bisection never fails.

## 2026-03-11 - CheckEncounter() Bug Fix

Properly ensure wrapping of GetDTforTrueAnomaly() return values onto [0,period].  This is clearly correct behavior and may fix bugs.

## 2026-03-08 - CheckEncounter() Bug Fixes

### `CheckEncounter` (PatchedConics.cs)

Both of the GetClosestApproach()+EncountersBody() codepaths are now protected by checks against the SOI radius so that if there's no chance of an intercept then the check will not run.  If the always-show-markers config is set and the celestial is the target then the markers will now always show up.  This fixes always-show-markers as well to always
do both checks and return the best one--and benefit from the initial guess fixes passed into GetClosestApproach().

## 2026-03-08 - CheckEncounter() bug Fixes

### `CheckEncounter` (PatchedConics.cs)

Compute the maxDT bounds passed into GetClosestApproach() for use in the Halley solver based on the orbital velocities at the geometric MOID point and the sphereOfInfluence, and seed the Halley solver with the UT of the geometric MOID point.  The previous values were some strange guesses based on the previously used bisection method.

## 2026-03-02 — Patched Conic Solver Bug Fixes

### `CheckEncounter` (PatchedConics.cs)

Restored the second-intercept fallback: if the first geometric intercept doesn't produce an SOI encounter, the solver is re-seeded near the second intercept and `EncountersBody` is called again. Also restored the `sqrt(SMA/SOI)` Pe/Ap buffer heuristic for targeted bodies. Kept try/catch, `logErrors` parameter, and thread-safety guards from the newer version.

### `GetClosestApproach` (PatchedConics.cs)

Reverted search window logic from synodic-period centering back to the original: full period for elliptical orbits, SOI-bounded for hyperbolic. Eliminates the dead code branch and the dependency on the previous `p.UTappr` value. Kept `pars.epsilon` from the newer API.

### `SolveClosestApproach` (Orbit.cs)

Removed the internal intercept-selection block that conflicted with `_CheckEncounter`'s two-call strategy. Fixed a sign error in the presolve quadratic formula (`drdv` → `-drdv`) and differentiated the two branches to select the appropriate root based on whether the bodies are approaching or receding.
