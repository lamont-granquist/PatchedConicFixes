# Stock intercept subdivision correction on PCF 0.1.4

## Scope and evidence

This candidate starts from official release `v0.1.4`, commit `135f3bb`, and carries forward the subdivision correction from `cdae23c` (the previously tested local 0.1.3.1 candidate). It is identified as local version `0.1.4.1`, informational version `0.1.4-intercept-candidate.1`, on branch `fix/intercept-subdivision-wrap-v014`.

The 10 October 2026 tourist capture trace shows `OrbitTargeter.UpdateNodesAndVectors -> Orbit._FindClosestPoints -> Targeting.Intercepts -> Targeting.add_crossing_subdivisions` on the main thread. At UT `28152013.389861476`, the same interval endpoints (`1.2999800827696877`, `0.099998467905360572`) recur at depths 64, 256, and 1024 while the interval count grows from 140 to 2060.

The stock `Targeting.Interval.Subdivide` unwraps a seam-crossing angular interval using a full turn, but wraps its two interior points using only a half turn. Its child intervals can expand and revisit the same range. The Harmony transpiler replaces precisely those two `pi` subtractions with `2*pi`. It validates both matches before making either change and retains KSP's sampling, pooling, linking, root search, and intersection return contract.

The PCF 0.1.4 encounter prefix, `HarmonyPatches.cs`, and `Baluev.cs` are unchanged from the official release. The new patch addresses stock targeting callers that still run outside PCF's Baluev encounter pipeline. The separate experimental cross-reference-body guard from local 0.1.3.2 is not included. Startup logs identify the candidate version and successful subdivision patch activation. Builds write to repository output only; deployment is disabled by default and no installation was performed.

## Offline validation completed on 11 October 2026

Validation uses the installed KSP 1.12.5 assemblies and a .NET Framework 4.8 test host.

- Unchanged official 0.1.4 baseline, with repository-only build configuration: 1289 tests, 1277 passed, 12 failed.
- Candidate: 1305 tests, 1293 passed, 12 failed. The failed test names match the baseline exactly; there are no additional failures. The upstream failures remain unresolved by this focused repair.
- All 16 targeting regression tests passed: angular partition/progress controls; the captured interval; a negative control proving that stock subdivision fails the progress assertion for that interval and the correction passes; analytic ellipse/circle intersections; saved Tourist A/B approach approximations; and the deterministic 200-case station geometry sweep.
- The full startup `PatchAll` integration test passed. It invokes the actual patched `PatchedConics._CheckEncounter`, prohibits fallback to `CheckEncounterOldMOID`, verifies the existing lunar encounter result, and also checks the stock targeting intersection result.
- Release build: zero warnings and zero errors. The separately built Release DLL and the DLL used by the tests have identical SHA-256 hashes.

Candidate DLL: `build/v014-release/PatchedConicFixes.dll`.

SHA-256: `5E9D597FC1172D6A56DB875E0490A65552BC73554517717AE445A346987E785E`.

Installed official 0.1.4 DLL remains unchanged, SHA-256: `1C3F1B86FBEB059B24A4500772D8F108846B5A10D985C4CDE5F5CDA002C494EA`.

TRX evidence is retained under `build/v014-validation/{baseline,targeting,candidate}`. The full captured orbit pair was not present in the hang record: the new regression reproduces the subdivision error with the exact interval endpoints, not the complete finite-burn hang. A repeat capture burn with the station target selected is still required to establish runtime success and useful next-orbit intercepts.

## Reproduction

From the repository root in PowerShell, with the existing local KSP path configuration:

```powershell
dotnet test Source/Tests/PatchedConicFixes.Tests.csproj -c Release -p:DeployToKSP=false "-p:OutputPath=$PWD/build/v014-candidate/" --filter FullyQualifiedName~TargetingSubdivisionTests
dotnet test Source/Tests/PatchedConicFixes.Tests.csproj -c Release --no-build -p:DeployToKSP=false "-p:OutputPath=$PWD/build/v014-candidate/"
dotnet build Source/PatchedConicFixes.csproj -c Release -p:DeployToKSP=false "-p:OutputPath=$PWD/build/v014-release/"
```

The full suite is expected to report the documented 12 baseline failures. Do not describe it as wholly passing.
