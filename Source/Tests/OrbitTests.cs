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

            var p         = new Orbit( 0.122026202733377, 0.971510969524243, 234243783.129783, 358.94830416185, 179.012674865338, 0.000732329708470401, 432322.256460584, earth);
            p.StartUT = 432322.256460584;
            p.EndUT   = 1560592.39540065;
            var nextPatch = new Orbit();
            var pars      = new PatchedConics.SolverParameters();

            bool result = HarmonyPatches.CheckEncounter(p, nextPatch,  432322.256460584, moon.orbitDriver, moon, pars, false);

            Assert.True(result);
        }

        [Fact]
        public void MoonTLIFailure1()
        {
            (CelestialBody earth, CelestialBody moon) = Bodies.MakeEarthMoon();

            var p         = new Orbit(0.122026203006779, 0.97151096952493, 234243783.135793, 358.948304186511, 179.012674841448, 0.000727985992839176, 432321.476460583, earth);
            var nextPatch = new Orbit();
            var pars      = new PatchedConics.SolverParameters();

            bool result = HarmonyPatches.CheckEncounter(p, nextPatch, 432321.476460583, moon.orbitDriver, moon, pars, false);

            Assert.True(result);
        }
    }
}
