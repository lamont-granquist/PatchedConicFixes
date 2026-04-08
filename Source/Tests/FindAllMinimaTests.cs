using Xunit;
using Xunit.Abstractions;

namespace PatchedConicFixes.Tests
{
    public class FindAllMinimaTests
    {
        private readonly ITestOutputHelper _testOutputHelper;

        public FindAllMinimaTests(ITestOutputHelper testOutputHelper)
        {
            _testOutputHelper = testOutputHelper;
        }

        [Fact]
        public unsafe void NearParabolicVesselTest()
        {
            Logger.Register(o => _testOutputHelper.WriteLine((string)o));

            // a = semi-major axis, e = eccentricity, i = inclination,
            // w = argument of periapsis, Om = longitude of ascending node
            var o1 = new COrbitData(386102172.159263, 0.99933864853071519, 1.6178598461249187, 6.2752580854970281, 0.31441741545137758);
            var o2 = new COrbitData(443719754.5067898, 0, 0, 0, 0);

            Baluev.MoidInfo* info = stackalloc Baluev.MoidInfo[4];

            int num = Baluev.FindAllMinima(o1, o2, info);

            num.ShouldEqual(3);
            info[0].dst.ShouldEqual(443464397.236679, 1e-10);
            info[0].u1.ShouldEqual(0.000288737023422094, 1e-10);
            info[0].u2.ShouldEqual(0.31404348382377, 1e-10);
            info[1].dst.ShouldEqual(10352230.1167095, 1e-10);
            info[1].u1.ShouldEqual(-1.72136528749143, 1e-10);
            info[1].u2.ShouldEqual(-2.82827362480136, 1e-10);
            info[2].dst.ShouldEqual(17379182.1969371, 1e-10);
            info[2].u1.ShouldEqual(1.72146689984871, 1e-10);
            info[2].u2.ShouldEqual(-2.82533036854555, 1e-10);
        }
    }
}
