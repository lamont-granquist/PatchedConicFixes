# Changelog

## 2026-03-02 — Patched Conic Solver Bug Fixes

### `_CheckEncounter` (PatchedConics.cs)

Restored the second-intercept fallback: if the first geometric intercept doesn't produce an SOI encounter, the solver is re-seeded near the second intercept and `EncountersBody` is called again. Also restored the `sqrt(SMA/SOI)` Pe/Ap buffer heuristic for targeted bodies. Kept try/catch, `logErrors` parameter, and thread-safety guards from the newer version.

### `_GetClosestApproach` (PatchedConics.cs)

Reverted search window logic from synodic-period centering back to the original: full period for elliptical orbits, SOI-bounded for hyperbolic. Eliminates the dead code branch and the dependency on the previous `p.UTappr` value. Kept `pars.epsilon` from the newer API.

### `_SolveClosestApproach` (Orbit.cs)

Removed the internal intercept-selection block that conflicted with `_CheckEncounter`'s two-call strategy. Fixed a sign error in the presolve quadratic formula (`drdv` → `-drdv`) and differentiated the two branches to select the appropriate root based on whether the bodies are approaching or receding.
