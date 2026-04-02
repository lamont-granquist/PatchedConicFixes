using Xunit;
using Xunit.Abstractions;

namespace PatchedConicFixes.Tests
{
    [Collection("KSP")]
    public class OrbitTests
    {
        private readonly ITestOutputHelper _testOutputHelper;

        public OrbitTests(ITestOutputHelper testOutputHelper)
        {
            _testOutputHelper = testOutputHelper;
        }

        // This is a bug with an typical TLI ejection to the Moon with the Moon ascending through a roughly equatorial
        // encounter, which exercises the need to test both MOID points and is buggy in KSP 1.12.x
        [Fact]
        public void AscendingMoonTLIFailure()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            (CelestialBody earth, CelestialBody moon) = Bodies.MakeEarthMoon();

            var p         = new Orbit(0.122026203095501, 0.971510969531635, 234243783.190063, 358.94830419701, 179.012674829019, 0.000723530899600804, 432320.676460583, earth) { StartUT = 432320.676460583, EndUT = 1560590.81583617 };
            var nextPatch = new Orbit();
            var pars      = new PatchedConics.SolverParameters();

            bool result = HarmonyPatches.CheckEncounter(p, nextPatch, 432320.676460583, moon.orbitDriver, moon, pars, false);

            Assert.True(result);

            Assert.Equal(1163413.6494660573, p.EndUT, 3);
            Assert.Equal(-6607296.1767470818, nextPatch.semiMajorAxis, 2);
            Assert.Equal(1.3779789903106174, nextPatch.eccentricity, 9);
            Assert.Equal(133.26553444141823, nextPatch.inclination, 7);
            Assert.Equal(205.38928286487098, nextPatch.LAN, 7);
            Assert.Equal(274.89492191030376, nextPatch.argumentOfPeriapsis, 7);
            Assert.Equal(-8.1599344523079314, nextPatch.meanAnomalyAtEpoch, 7);
            Assert.Equal(1163413.6494660573, nextPatch.epoch, 3);

            Vector3d vesselPos = p.getPositionAtUT(p.EndUT);
            Vector3d moonPos   = moon.orbitDriver.orbit.getPositionAtUT(p.EndUT);

            Assert.Equal(moon.sphereOfInfluence, (vesselPos - moonPos).magnitude, 1);
        }

        // This is the same case as AscendingMoonTLIFailure() but we've warped past the first MOID point and
        // even with the buggy solver the encounter will appear.
        [Fact]
        public void AscendingMoonTLISuccessful()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            (CelestialBody earth, CelestialBody moon) = Bodies.MakeEarthMoon();

            var p         = new Orbit(0.122026202811084, 0.97151096953721, 234243783.237659, 358.948304503794, 179.012674522348, 0.00924041893390898, 433850.052085614, earth) { StartUT = 433850.052085614, EndUT = 1562120.19180508 };
            var nextPatch = new Orbit();
            var pars      = new PatchedConics.SolverParameters { TimeSolverIterations1 = 3 };

            bool result = HarmonyPatches.CheckEncounter(p, nextPatch, 433850.052085614, moon.orbitDriver, moon, pars, false);

            Assert.True(result);
            Assert.Equal(1163413.6495317207, p.EndUT, 3);
            Assert.Equal(-6607296.1804669211, nextPatch.semiMajorAxis, 2);
            Assert.Equal(1.3779790095658691, nextPatch.eccentricity, 9);
            Assert.Equal(133.26553507031255, nextPatch.inclination, 7);
            Assert.Equal(205.38928212367242, nextPatch.LAN, 7);
            Assert.Equal(274.89492052411111, nextPatch.argumentOfPeriapsis, 7);
            Assert.Equal(-8.1599344587965881, nextPatch.meanAnomalyAtEpoch, 6);
            Assert.Equal(1163413.6495317207, nextPatch.epoch, 3);

            Vector3d vesselPos = p.getPositionAtUT(p.EndUT);
            Vector3d moonPos   = moon.orbitDriver.orbit.getPositionAtUT(p.EndUT);

            Assert.Equal(moon.sphereOfInfluence, (vesselPos - moonPos).magnitude, 0);
        }

        // This is the same case as AscendingMoonTLIFailure() but we've warped past the first MOID point and
        // even with the buggy solver the encounter will appear.  This exercises a bug in the SolveSOI() Newton
        // solver where it runs into a "flat" region and the denominator goes to zero and the step becomes ~1.8e+19
        [Fact]
        public void FlickeringTyloEncounter()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            (CelestialBody jool, CelestialBody tylo) = Bodies.MakeJoolTylo();

            var p         = new Orbit(0.841141131311559, 0.0263025544044183, 64166196.6687604, 263.186581736814, 131.032510427265, 1.63048082893379, 1261109.62120736, jool) { StartUT = 1261109.62120736, EndUT = 1453245.5413889 };
            var nextPatch = new Orbit();
            var pars      = new PatchedConics.SolverParameters { TimeSolverIterations1 = 3 };

            bool result = HarmonyPatches.CheckEncounter(p, nextPatch, 1261109.62120736, tylo.orbitDriver, tylo, pars, false);

            Assert.True(result);
            Assert.Equal(1409549.6282595578, p.EndUT, 3);
            Assert.Equal(-3070988.0775043946, nextPatch.semiMajorAxis, 1);
            Assert.Equal(2.7065531829057434, nextPatch.eccentricity, 7);
            Assert.Equal(4.8387567394518829, nextPatch.inclination, 6);
            Assert.Equal(141.84225399409706, nextPatch.LAN, 6);
            Assert.Equal(218.03375183903731, nextPatch.argumentOfPeriapsis, 5);
            Assert.Equal(-2.5337155636348196, nextPatch.meanAnomalyAtEpoch, 7);
            Assert.Equal(1409549.6282595578, nextPatch.epoch, 3);

            Vector3d vesselPos = p.getPositionAtUT(p.EndUT);
            Vector3d moonPos   = tylo.orbitDriver.orbit.getPositionAtUT(p.EndUT);

            Assert.Equal(tylo.sphereOfInfluence, (vesselPos - moonPos).magnitude, 1);
        }
    }
}
