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

        [Fact]
        public void NearParabolicVesselTest()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            // a = semi-major axis, e = eccentricity, i = inclination,
            // w = argument of periapsis, Om = longitude of ascending node
            var o1 = new COrbitData(386102172.159263, 0.99933864853071519, 1.6178598461249187, 6.2752580854970281, 0.31441741545137758);
            var o2 = new COrbitData(443719754.5067898, 0, 0, 0, 0);

            SMOIDResult result = Baluev.MOID_fast(in o1, in o2);

            //Assert.True(result.good, "Result should be marked as reliable");

            result.distance.ShouldEqual( 10352230.116709549, 1e-14);
            result.u1.ShouldEqual(-1.7213652874914271, 1e-15);
            result.u2.ShouldEqual( -2.8282736248013554, 1e-15);
            result.rootCount.ShouldEqual(12);
            result.distanceError.ShouldEqual(2.787539797647343e-07, 1e-15);
            result.u1Error.ShouldEqual(1.0876683191309692e-15, 1e-15);
            result.u2Error.ShouldEqual(9.783818447810477e-16, 1e-15);
        }
    }
}
