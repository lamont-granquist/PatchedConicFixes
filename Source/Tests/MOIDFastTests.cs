using Xunit;
using Xunit.Abstractions;

namespace PatchedConicFixes.Tests
{
    public class MOIDFastTests
    {
        private readonly ITestOutputHelper _testOutputHelper;

        public MOIDFastTests(ITestOutputHelper testOutputHelper)
        {
            _testOutputHelper = testOutputHelper;
        }

        [Fact]
        public void EarthLikeAndEccentricOrbit()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            // a = semi-major axis, e = eccentricity, i = inclination,
            // w = argument of periapsis, Om = longitude of ascending node
            var o2 = new COrbitData(482701874.555625, 0.0, 0.0, 0.0, 0.0);
            var o1 = new COrbitData(466885815.754593, 0.845035162370339,
                1.13266029679883, 3.70233801890817, 1.22862857305205);

            SMOIDResult result = Baluev.MOID_fast(in o1, in o2);

            Assert.True(result.good, "Result should be marked as reliable");

            result.distance.ShouldEqual(6675627.898447883, 1e-14);
            result.u1.ShouldEqual(1.6023861460952622, 1e-15);
            result.u2.ShouldEqual(1.2342501326931417, 1e-15);
            result.rootCount.ShouldEqual(8);
            result.distanceError.ShouldEqual(3.021029828744109e-07, 1e-15);
            result.u1Error.ShouldEqual(9.910191270109278e-16, 1e-15);
            result.u2Error.ShouldEqual(1.0581479243090316e-15, 1e-15);
        }
    }
}
