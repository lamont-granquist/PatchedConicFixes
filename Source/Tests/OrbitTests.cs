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

        [Fact]
        public void MoonTLISuccessful()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            (CelestialBody earth, CelestialBody moon) = Bodies.MakeEarthMoon();

            var p         = new Orbit(0.122026202811084, 0.97151096953721, 234243783.237659, 358.948304503794, 179.012674522348, 0.00924041893390898, 433850.052085614, earth) { StartUT = 433850.052085614, EndUT = 1562120.19180508 };
            var nextPatch = new Orbit();
            var pars      = new PatchedConics.SolverParameters { TimeSolverIterations1 = 3 };

            bool result = HarmonyPatches.CheckEncounter(p, nextPatch, 433850.052085614, moon.orbitDriver, moon, pars, false);

            Assert.True(result);
            Assert.Equal(1163413.6495317207, p.EndUT);
            Assert.Equal(-6607296.1804669211, nextPatch.semiMajorAxis);
            Assert.Equal(1.3779790095658691, nextPatch.eccentricity);
            Assert.Equal(133.26553507031255, nextPatch.inclination);
            Assert.Equal(205.38928212367242, nextPatch.LAN);
            Assert.Equal(274.89492052411111, nextPatch.argumentOfPeriapsis);
            Assert.Equal(-8.1599344587965881, nextPatch.meanAnomalyAtEpoch);
            Assert.Equal(1163413.6495317207, nextPatch.epoch);

            Vector3d vesselPos = p.getPositionAtUT(p.EndUT);
            Vector3d moonPos   = moon.orbitDriver.orbit.getPositionAtUT(p.EndUT);

            Assert.Equal(moon.sphereOfInfluence, (vesselPos - moonPos).magnitude, 3);
        }

        [Fact]
        public void MoonTLIFailure()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            (CelestialBody earth, CelestialBody moon) = Bodies.MakeEarthMoon();

            var p         = new Orbit(0.122026203095501, 0.971510969531635, 234243783.190063, 358.94830419701, 179.012674829019, 0.000723530899600804, 432320.676460583 , earth) { StartUT = 432320.676460583, EndUT = 1560590.81583617 };
            var nextPatch = new Orbit();
            var pars      = new PatchedConics.SolverParameters();

            bool result = HarmonyPatches.CheckEncounter(p, nextPatch, 432320.676460583, moon.orbitDriver, moon, pars, false);

            Assert.True(result);

            Assert.Equal(1163413.6494660573, p.EndUT);
            Assert.Equal(-6607296.1767470818, nextPatch.semiMajorAxis);
            Assert.Equal(1.3779789903106174, nextPatch.eccentricity);
            Assert.Equal(133.26553444141823, nextPatch.inclination);
            Assert.Equal(205.38928286487098, nextPatch.LAN);
            Assert.Equal(274.89492191030376, nextPatch.argumentOfPeriapsis);
            Assert.Equal(-8.1599344523079314, nextPatch.meanAnomalyAtEpoch);
            Assert.Equal(1163413.6494660573, nextPatch.epoch);

            Vector3d vesselPos = p.getPositionAtUT(p.EndUT);
            Vector3d moonPos   = moon.orbitDriver.orbit.getPositionAtUT(p.EndUT);

            Assert.Equal(moon.sphereOfInfluence, (vesselPos - moonPos).magnitude, 3);
        }
    }
}
