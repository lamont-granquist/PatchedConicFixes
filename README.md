# PatchedConicFixes

This is a near overhaul of KSP 1.12.x's encounter resolution in the patched conic
solver.  It fixes several regressions introduced around KSP 1.8.0, along with
multiple design problems that KSP had from the start.  Many or most of the bugs
where SOI encounters to go undetected, even when your ship is on a clear intercept
course, have been fixed.  A few examples of affected scenarios include problems
with flickering orbits around Jool, and transfers to an ascending Moon in
RealSolarSystem, and many others. The replacement encounter solver is also
approximately 5× faster than the stock implementation.

## Known Limitations

This implementation has known failure conditions that stem from fundamental design
choices in the algorithm (shared with the old KSP implementation).  A meaningful fix
would require a ground-up rewrite (which is just dependent upon enough free time).

If you encounter still buggy behavior, it's quite likely a known issue. For the time being
bug reports for this mod are unnecessary.

## Requirements

- Kerbal Space Program 1.12.x
- KSPBurst (1.7.2 or compatible)
- Harmony (2.2.1 or compatible)

## Installation

Install via CKAN, or manually copy `GameData/PatchedConicFixes` into your KSP `GameData` folder.

## Building from Source

### Prerequisites

- [.NET SDK 6.0+](https://dotnet.microsoft.com/download)
- KSP 1.12.x installed

### Setup

Edit `Source/Directory.Build.props` to point to your KSP installation:

```xml
<KSPData>/path/to/KSP.app/Contents/Resources/Data</KSPData>     <!-- Mac -->
<KSPGameData>/path/to/KSP/GameData</KSPGameData>
```

Or create a `Source/Directory.Build.props.user` (gitignored) to keep your local paths out of version control:

```xml
<Project>
  <PropertyGroup>
    <KSPData>/your/local/path/KSP.app/Contents/Resources/Data</KSPData>
    <KSPGameData>/your/local/path/GameData</KSPGameData>
  </PropertyGroup>
</Project>
```

### Build

```
cd Source
dotnet build
```

The compiled DLL will be placed in `GameData/PatchedConicFixes/Plugins/`.

## License

MIT

## Changes to KSP APIs

### Orbit.FindClosestPoints() and the Targetting class

This algorithm finds the closest points between two conic sections.  It was entirely replaced by an algorithm
by Baluev and Mikryukov:

- Baluev RV, Mikryukov DV. Fast error-controlling MOID computation for confocal elliptic orbits. Astronomy and Computing. 2019 Apr;27:11–22. doi:10.1016/j.ascom.2019.02.005
- Mikryukov D, Baluev R. A lower bound of the distance between two elliptic orbits. Celest Mech Dyn Astr. 2019 Jun;131(6):28. doi:10.1007/s10569-019-9907-3
- Baluev RV. Fast error-safe MOID computation involving hyperbolic orbits. Astronomy and Computing. 2021 Jan;34:100440. doi:10.1016/j.ascom.2020.100440

The c++ code from distlink.cpp was converted to burst-compiled C# code.  It constructs the 16th degree complex polynomial via the DFT in those
papers and solves it via Newton root finding and deflation using the tuned initial guesses from those papers.  It has been modified to return up to
6 minima points unlike the stock routine which only returns 2 points and can miss some minima due to its uniform sampling approach.

I believe there can be only 4 actual minima, but up to 6 "distinct" points were returned in some cases, likely due to degenerate roots which escaped
the de-duplication checker.  I don't think more than 12 real roots can exist in this problem so 6 should be a correct upper limit.

Since we definitely need at least 3 minima, even with an SOI filter, the old FindClosestPoints() API couldn't be fixed.  The code is still there, it
still behaves the same way, there are two other references to it in the codebase.

### CheckEncounter()

This now runs through all minima returned by the Baluev algorithm that survive a filter on the SOI radius, ordered by transit time to the minima point
from the current vessel position.  The alogorithm additional tries the last one first (minus a hole period) and the first one last (plus a whole period).
This can find some crossings where the MOID point is e.g. back in time, but the time-domain actual SOI crossing still lies in the future of the Vessel.
Doing this work is necessary for correctness, every possible geometrical approach should be checked in order.

### GetClosestApproach()

The validation in this API has been removed.  It validated that the geometrical MOID points were at valid times on the patch (e.g. up to one period in the
future on elliptical orbits).  That validation has be moved to after the SOI crossing is actually found.  That work needs to get done due to correctness,
even though it might get thrown away by the lazy instead of eager validation check.

### SolveClosestApproach()

This is the routine which takes a geometrical closest approach and finds the actual time-domain closest approach between the Vessel and the CelestialBody on
the actual Orbits.  It has been changed to do an expanding search around the MOID point to find a sign-change in the derivative.  It then does a bracketed
Halley search with bisection as a fallback to avoid cycles and other pathologies.  A lot of the presolving steps have been removed from this algorithm in
favor of the initial bracketing.

### EncountersBody()

This is the routine which drives taking the closest approach, calling SolveSOI() and turning it into an SOI crossing and updating the actual Orbit patch.
This now calls the new internal ValidateSOICrossingTime() which does the work that GetClosestApproach() used to do on the MOID time, but applies it to the
solved SOI time instead.

### SolveSOI()

This routine has been changed to be a bracketing Newton search with fallback to bisection to remove cycles and other pathologies common in Newton's method,
the metric was also changed to the square of the distance miss (which should behave better far from the solution) and it calls a cheaper routine to calculate
the derivative function.  This was validated to fix some real world reported orbit-flickering issues.

## Future Work

The whole MOID-solving approach probably should be rejected, or kept only as a boolean prefilter.

The biggest problem is orbits where the Vessel is nearly in the same Orbit as the CelestialBody.  In that case you could have every point in the Vessel Orbit
potentially close enough to be within the SOI of the CelestialBody, only with slightly different periods.  In that case, you could get an SOI crossing anywhere
and the MOID gives no information and seeding guesses at the MOID points is likely to miss actual time-domain SOI crossings.

What should probably happen instead is implementing the cross track and synodic prefilters from:

- Ogle A, Spitas C. A Novel Analytical Method to Determine Future Close Approaches between Satellites. 2022.

Then working directly on the time-domain distance function between the actual two orbits, approximate them with Chebyshev polynomials and solve for where
the derivative crosses zero (replacing the algorithm in SolveClosestApproach()) using the algorithm in:

- Denenberg E. Satellite closest approach calculation through Chebyshev Proxy Polynomials. Acta Astronautica. 2020 May;170:55–65. doi:10.1016/j.actaastro.2020.01.020

