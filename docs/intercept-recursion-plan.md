# Intercept recursion investigation and fix plan

## Objective

Make the intersection/closest-points calculation terminate correctly for the observed Mun encounter/burn scenario, while retaining useful intersection display and PatchedConicFixes' encounter detection.

## Established baseline

- Fork: `maddavo/PatchedConicFixes`; upstream: `lamont-granquist/PatchedConicFixes`.
- Investigation branch: `fix/intercept-recursion`, starting at `c6642ac`.
- Installed DLL and repository DLL are identical, assembly version `0.1.3.0`.
- SHA-256: `07E8FF25A297952CF84ACAA9388879911A337F63C3EA4736E44D50CBE62967D8`.
- The current `_CheckEncounter` prefix calls `CheckEncounterOldMOID`; the newer `CheckEncounter` invocation is commented out. Do not silently switch solvers as part of this investigation.
- No repository AGENTS.md was found during initial inspection.

## Evidence and limits

On 9 October 2026, the caller trace recorded the active chain:

`OrbitUtil.CalculateIntersections -> Orbit._FindClosestPoints -> Targeting.Intercepts`

`Targeting.add_crossing_subdivisions` reached depths 64, 256, and 1024, with output interval counts 139, 523, and 2059. The scene was FLIGHT and the active vessel was Tourist A. The player was executing a burn with intersections visible, alternating manual control and MechJeb.

An earlier failure explicitly logged a StackOverflowException with repeated `Targeting.Interval.FindRoot` frames. Treat both recursive paths as investigation targets; do not claim they have the same numerical trigger until reproduced.

The trace identifies the active intersection route, but does not contain the exact input orbit pair, interval anomaly endpoints, or the outer ManeuverNodeEditor method. A nearby PCF encounter-check message is not proof that its call caused this recursive event.

Upstream issue: https://github.com/lamont-granquist/PatchedConicFixes/issues/2

The shipped KSP method inspection shows crossing subdivision recursively follows an interval whose endpoint inside/outside classifications disagree. That method has no explicit depth or interval-width termination guard. This identifies a missing safeguard; the reason the observed interval did not converge still requires numerical evidence.

## Implementation sequence

1. **Make builds safe and reproducible.** Change deployment to an explicit opt-in operation before building: the existing project automatically copies every build into the configured KSP GameData. Configure local KSP_x64_Data paths and correct the test project's hard-coded Burst/Mathematics paths. Record compiler and dependency versions.
2. **Create a replay case.** Recover orbit inputs from the authoritative save and saved evidence where possible. Capture the exact reference/target orbit pair and anomaly endpoints only if those cannot be recovered. Investigate angular wrapping, invalid samples, and whether a selected child interval actually shrinks. Do not request another large dump merely to establish the already-recorded caller path.
3. **Implement the numerical fix.** Patch the narrow stock routine through Harmony within this fork. Require finite inputs, representable progress, and bounded work. Address both crossing subdivision and the invalid-sample recursive branch in FindRoot where warranted. Preserve sample/interval pool ownership and linked interval ordering. A depth cap alone is containment, not proof of corrected intersections. Define explicit failure behavior so unresolved calculations cannot become fabricated valid intersections.
4. **Verify calculations offline.** Add a regression for the captured failure and representative elliptical/hyperbolic, wrap-boundary, near-tangent, and invalid-region cases. Check termination, finite results, intersection accuracy against an independent reference, and normal SOI encounter behavior. Run the existing relevant solver tests. Build/test success is not flight validation.
5. **Validate in KSP.** Back up the installed DLL, install a uniquely identified candidate while KSP is closed, and verify its hash. Use a copy of the authoritative save to repeat the target/intersection/burn scenario. Confirm useful intersections and encounter prediction remain correct, not merely that the game survives. Retain focused diagnostic logging for this validation.
6. **Prepare the contribution.** After runtime validation, prepare a focused upstream PR referring to issue #2, documenting the trigger, fix, numerical tests, and flight results. Keep proprietary KSP DLLs, personal saves, full dumps, and machine paths out of the contribution. Posting an upstream comment/report requires Dave's explicit authorization.

## Completion criteria

- The captured numerical case terminates and produces validated results or an explicit unresolved result.
- Ordinary intersections and SOI encounter detection remain accurate in regression checks.
- The original gameplay sequence succeeds with the candidate loaded and verified.
- The change and rollback artifact are reviewable and ready for an upstream contribution.

## Investigation progress

- Quicksave #70 is the authoritative starting state. A read-only replay copy is retained in the ignored build directory. SHA-256: `0AE5848563C581EF4E308892971BD5EB8DBE0CDAD45B4330E6F4FECF6EDF36B7`.
- Both tourists are on Mun approach trajectories. Tourist A targets Mun Alpha Station and has a saved periapsis node at UT `27642069.704133388`, with approximately `266.81309009934671 m/s` retrograde delta-v. Tourist B has no saved maneuver node in this save.
- The stock `Targeting.Interval.Subdivide` subtracts pi instead of two-pi from each interior true anomaly when crossing the +pi/-pi seam. Real KSP method tests reproduce a child width of `3.1612300408133764` radians from a parent only `0.083185307179586054` radians wide. A correction of the two wrap subtraction constants restores shrinking, ordered children and preserves total interval width.
- Seam and non-seam tests, an independent analytic ellipse/circle crossing check, Tourist A impulse-burn approximations, Tourist B's starting orbit, and a 200-case same-body station sweep exercise the candidate. The finite-burn hang itself has not been reproduced offline: no trace record yet contains the precise hung orbit pair and anomaly endpoints.
- The Intercept UI calls CalculateIntersections while walking escape/encounter patches before checking whether the patch and target have the same reference body. Returning default/no intersections for incompatible frames is a separate candidate change, so its effect can be measured independently from the wrap correction.
- Unmodified upstream tests on this machine: 1289 total, 1277 passing, 12 failing. The first combined candidate run had exactly the same 12 failure names and 14 additional passing regression cases. Ten failures are explicitly named BrokenTestsTodo; the other two are FindAllMinimaTests.NearParabolicVesselTest and SqDistSurfaceTests.SqDist2DIter_StepSizeDecreasesNearMinimum. Do not describe the full suite as passing.
- Test recursion/sample budgets are sentries to protect the test host. They do not ship in the mod. The numerical candidate corrects interval wrapping; it does not establish a general bound on every remaining stock search path.
- Keep separate commits/candidates for the wrap correction and cross-frame guard. Start runtime validation with the wrap correction alone and the enhanced trace. Completion still requires useful next-orbit station intercepts and a successful finite capture burn, plus investigation of any repeated failure with exact captured inputs.
