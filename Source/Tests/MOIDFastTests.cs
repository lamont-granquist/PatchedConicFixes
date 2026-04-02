using System;
using System.Diagnostics;
using Xunit;
using Xunit.Abstractions;

namespace PatchedConicFixes.Tests
{
    public unsafe class MOIDFastTests
    {
        private readonly ITestOutputHelper _testOutputHelper;

        public MOIDFastTests(ITestOutputHelper testOutputHelper)
        {
            _testOutputHelper = testOutputHelper;
        }

        static void AssertClose(double expected, double actual, double tol, string label = "")
        {
            Assert.True(Math.Abs(expected - actual) < tol,
                $"{label} Expected {expected:G17}, got {actual:G17}, diff = {Math.Abs(expected - actual):G17}");
        }

        [Fact]
        public void EarthLikeAndEccentricOrbit()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            // a = semi-major axis, e = eccentricity, i = inclination,
            // w = argument of periapsis, Om = longitude of ascending node
            var O2 = new COrbitData(482701874.555625, 0.0, 0.0, 0.0, 0.0);
            var O1 = new COrbitData(466885815.754593, 0.845035162370339,
                1.13266029679883, 3.70233801890817, 1.22862857305205);

            SMOIDResult result = Baluev.MOID_fast(in O1, in O2);

            Assert.True(result.good, "Result should be marked as reliable");

            AssertClose(6675627.898447883, result.distance, 1e-7, "distance");
            AssertClose(1.6023861460952622, result.u1, 1e-15, "u1");
            AssertClose(1.2342501326931417, result.u2, 1e-15, "u2");
            Assert.Equal(8, result.rootCount);
            AssertClose(3.021029828744109e-07, result.distanceError,  1e-15, "distanceError");
            AssertClose(9.910191270109278e-16, result.u1Error, 1e-15, "u1Error");
            AssertClose(1.0581479243090316e-15, result.u2Error, 1e-15, "u2Error");
        }
    }
}
