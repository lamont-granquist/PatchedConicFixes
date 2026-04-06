using System;
using Xunit;
using Xunit.Abstractions;

namespace PatchedConicFixes.Tests
{
    public unsafe class NewtonSqDistTests
    {
        private readonly ITestOutputHelper _testOutputHelper;

        public NewtonSqDistTests(ITestOutputHelper testOutputHelper)
        {
            _testOutputHelper = testOutputHelper;
        }

        private const double Tol = 1e-10;

        private static void AssertClose(double expected, double actual, double tol = Tol, string label = "")
        {
            Assert.True(Math.Abs(expected - actual) < tol,
                $"{label} Expected {expected:G17}, got {actual:G17}, diff = {Math.Abs(expected - actual):G17}");
        }

        // --- basic convergence ---

        [Fact]
        public void ConvergesToCriticalPoint()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 1.0, u2 = 0.5;
            Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out double u1Err, out double u2Err,
                out double rho, out double rhoErr,
                1e-14, out ulong iterCnt, 30, g, H);

            double gNorm = Math.Sqrt(g[0] * g[0] + g[1] * g[1]);
            Assert.True(gNorm < 1e-10,
                $"Gradient should be near zero, got {gNorm:G17}");
        }

        [Fact]
        public void ReturnsPositiveDefiniteAtMinimum()
        {
            // Start near a minimum, should converge and report +2
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            // Use MinDistanceFor to get a good starting point near a minimum
            double u1Start = 1.0;
            Baluev.MinDistanceFor(&data, u1Start, out double u2Start);
            double u1 = u1Start, u2 = u2Start;

            short hsign = Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out _, out _,
                1e-14, out _, 30, g, H);

            // Should be a minimum (+2) or at least positive semi-definite (+1)
            Assert.True(hsign >= 1,
                $"Expected positive (semi-)definite Hessian, got {hsign}");
        }

        [Fact]
        public void RhoConsistentWithDistanceBetween()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 1.0, u2 = 0.5;
            Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out double rho, out _,
                1e-14, out _, 30, g, H);

            double dist        = Baluev.DistanceBetween(&data, true, true, u1, u2);
            double distFromRho = Math.Sqrt(2.0 * data.a1 * data.a2 * rho);
            AssertClose(dist, distFromRho, 1e-8, "distance from rho");
        }

        // --- error estimates ---

        [Fact]
        public void ErrorEstimatesAreNonNegative()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 1.0, u2 = 0.5;
            Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out double u1Err, out double u2Err,
                out _, out double rhoErr,
                1e-14, out _, 30, g, H);

            Assert.True(u1Err >= 0, $"u1Err should be ≥ 0, got {u1Err:G17}");
            Assert.True(u2Err >= 0, $"u2Err should be ≥ 0, got {u2Err:G17}");
            Assert.True(rhoErr >= 0, $"rhoErr should be ≥ 0, got {rhoErr:G17}");
        }

        [Fact]
        public void ErrorEstimatesAreSmallAtConvergence()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1Start = 1.0;
            Baluev.MinDistanceFor(&data, u1Start, out double u2Start);
            double u1 = u1Start, u2 = u2Start;

            Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out double u1Err, out double u2Err,
                out _, out double rhoErr,
                1e-14, out _, 30, g, H);

            Assert.True(u1Err < 1e-8,
                $"u1Err should be small at convergence, got {u1Err:G17}");
            Assert.True(u2Err < 1e-8,
                $"u2Err should be small at convergence, got {u2Err:G17}");
            Assert.True(rhoErr < 1e-10,
                $"rhoErr should be small at convergence, got {rhoErr:G17}");
        }

        // --- iteration count ---

        [Fact]
        public void IterationCountIsReasonable()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 1.0, u2 = 0.5;
            Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out _, out _,
                1e-14, out ulong iterCnt, 30, g, H);

            Assert.True(iterCnt <= 30,
                $"Should converge within maxcount, got {iterCnt}");
        }

        [Fact]
        public void AlreadyAtMinimum_FewIterations()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            // Pre-converge with several SqDist2DIter steps
            double u1 = 1.0, u2 = 0.5;
            for (int i = 0; i < 20; i++)
                Baluev.SqDist2DIter(&data, true, true, ref u1, ref u2,
                    out _, g, H, out _);

            double u1Pre = u1, u2Pre = u2;
            Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out _, out _,
                1e-14, out ulong iterCnt, 30, g, H);

            Assert.True(iterCnt <= 5,
                $"Starting near minimum should converge quickly, got {iterCnt}");
            AssertClose(u1Pre, u1, 1e-10, "u1 shouldn't move much");
            AssertClose(u2Pre, u2, 1e-10, "u2 shouldn't move much");
        }

        // --- angle wrapping ---

        [Fact]
        public void AngleWrapping_ResultInRange()
        {
            // Start with u1 near ±π boundary
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 3.0, u2 = -3.0; // near ±π
            Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out _, out _,
                1e-14, out _, 30, g, H);

            Assert.True(u1 >= -Math.PI - 0.01 && u1 <= Math.PI + 0.01,
                $"u1 = {u1:G17} should be near [-π, π]");
            Assert.True(u2 >= -Math.PI - 0.01 && u2 <= Math.PI + 0.01,
                $"u2 = {u2:G17} should be near [-π, π]");
        }

        // --- different orbit configurations ---

        [Fact]
        public void CoplanarCircularOrbits()
        {
            var     O1   = new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0);
            var     O2   = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 0.5, u2 = 0.5;
            short hsign = Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out double rho, out _,
                1e-14, out _, 30, g, H);

            // Minimum distance between concentric circles is |a2-a1| = 1
            double dist = Math.Sqrt(2.0 * data.a1 * data.a2 * rho);
            AssertClose(1.0, dist, 1e-6, "concentric circle distance");

            // u1 should equal u2 at minimum (same angle)
            AssertClose(u1, u2, 1e-6, "u1 = u2 at minimum");
        }

        [Fact]
        public void PerpendicularPlanes()
        {
            var     O1   = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var     O2   = new COrbitData(2.0, 0.0, Math.PI / 2.0, 0.0, 0.0);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            // Start near the intersection line
            double u1 = 0.1, u2 = 0.1;
            short hsign = Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out double rho, out _,
                1e-14, out _, 30, g, H);

            // Same radius circles with perpendicular planes intersect: MOID = 0
            double dist = Math.Sqrt(2.0 * data.a1 * data.a2 * rho);
            AssertClose(0.0, dist, 1e-6, "intersecting orbits MOID");
        }

        [Fact]
        public void HighEccentricityOrbits()
        {
            var     O1   = new COrbitData(2.0, 0.8, 0.3, 0.5, 0.0);
            var     O2   = new COrbitData(3.0, 0.6, 0.6, 1.0, 0.8);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 0.5;
            Baluev.MinDistanceFor(&data, u1, out double u2);

            short hsign = Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out double rho, out _,
                1e-14, out _, 30, g, H);

            double gNorm = Math.Sqrt(g[0] * g[0] + g[1] * g[1]);
            Assert.True(gNorm < 1e-8, $"Gradient should vanish, got {gNorm:G17}");
            Assert.True(rho >= 0, $"rho should be non-negative, got {rho:G17}");
        }

        // --- multiple starting points converge to critical points ---

        [Fact]
        public void MultipleStartingPoints_FiniteResults()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            var starts = new (double u1, double u2)[] { (0.5, 0.5), (1.0, -1.0), (2.0, 0.0), (-1.5, 1.5), (0.0, 2.5), (-2.0, -2.0) };

            foreach ((double u1Start, double u2Start) in starts)
            {
                double u1 = u1Start, u2 = u2Start;
                Baluev.NewtonSqDist(&data, ref u1, ref u2,
                    out _, out _, out double rho, out _,
                    1e-14, out _, 30, g, H);

                Assert.True(!double.IsNaN(rho) && !double.IsInfinity(rho),
                    $"Starting ({u1Start:F1},{u2Start:F1}): rho = {rho:G17}");
                Assert.True(!double.IsNaN(u1) && !double.IsNaN(u2),
                    $"Starting ({u1Start:F1},{u2Start:F1}): u1={u1:G17}, u2={u2:G17}");
            }
        }

        [Fact]
        public void SeededStartingPoints_AllConverge()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            for (int k = 0; k < 8; k++)
            {
                double u1 = -Math.PI + 0.1 + k * Math.PI / 4.0;
                double d  = Baluev.MinDistanceFor(&data, u1, out double u2);
                if (d < 0) continue;

                Baluev.NewtonSqDist(&data, ref u1, ref u2,
                    out _, out _, out double rho, out _,
                    1e-14, out _, 30, g, H);

                double gNorm = Math.Sqrt(g[0] * g[0] + g[1] * g[1]);
                _testOutputHelper.WriteLine(gNorm.ToString());
                Assert.True(gNorm < 0.1,
                    $"Seeded from u1={-Math.PI + 0.1 + k * Math.PI / 4.0:F3}: gradient = {gNorm:G6}");
            }
        }

        // --- Hessian sign classification ---

        [Fact]
        public void HessianSign_IsValidValue()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            var starts = new (double u1, double u2)[] { (0.5, 0.5), (1.0, -1.0), (2.0, 0.0) };

            foreach ((double u1Start, double u2Start) in starts)
            {
                double u1 = u1Start, u2 = u2Start;
                short hsign = Baluev.NewtonSqDist(&data, ref u1, ref u2,
                    out _, out _, out _, out _,
                    1e-14, out _, 30, g, H);

                Assert.True(hsign >= -2 && hsign <= 2,
                    $"Hsign should be in [-2,+2], got {hsign}");
            }
        }

        [Fact]
        public void HessianSign_PositiveDefiniteConsistentWithMinimum()
        {
            // When hsign = +2, rho should be a local minimum:
            // perturbing u1 or u2 should increase rho
            var     O1     = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2     = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data   = SAuxData.Create(in O1, in O2);
            double* g      = stackalloc double[2];
            double* H      = stackalloc double[3];
            double* gDummy = stackalloc double[2];

            double u1 = 1.0;
            Baluev.MinDistanceFor(&data, u1, out double u2);

            short hsign = Baluev.NewtonSqDist(&data, ref u1, ref u2,
                out _, out _, out double rho, out _,
                1e-14, out _, 30, g, H);

            if (hsign != +2) return; // only test if we actually found a minimum

            double h     = 1e-5;
            double rhoP1 = Baluev.SqDistVG(&data, true, true, u1 + h, u2, gDummy);
            double rhoM1 = Baluev.SqDistVG(&data, true, true, u1 - h, u2, gDummy);
            double rhoP2 = Baluev.SqDistVG(&data, true, true, u1, u2 + h, gDummy);
            double rhoM2 = Baluev.SqDistVG(&data, true, true, u1, u2 - h, gDummy);

            Assert.True(rhoP1 >= rho - 1e-12, $"rho(u1+h) = {rhoP1:G17} < rho = {rho:G17}");
            Assert.True(rhoM1 >= rho - 1e-12, $"rho(u1-h) = {rhoM1:G17} < rho = {rho:G17}");
            Assert.True(rhoP2 >= rho - 1e-12, $"rho(u2+h) = {rhoP2:G17} < rho = {rho:G17}");
            Assert.True(rhoM2 >= rho - 1e-12, $"rho(u2-h) = {rhoM2:G17} < rho = {rho:G17}");
        }

        // --- symmetry ---

        [Fact]
        public void SymmetricOrbits_SymmetricResult()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.0, 0.0);
            var     O2   = new COrbitData(2.0, 0.3, 0.5, Math.PI, 0.0);
            var     data = SAuxData.Create(in O1, in O2);
            double* g1   = stackalloc double[2];
            double* H1   = stackalloc double[3];
            double* g2   = stackalloc double[2];
            double* H2   = stackalloc double[3];

            double u1a = 0.5, u2a = 0.5;
            Baluev.NewtonSqDist(&data, ref u1a, ref u2a,
                out _, out _, out double rhoA, out _,
                1e-14, out _, 30, g1, H1);

            double u1b = -0.5, u2b = -0.5;
            Baluev.NewtonSqDist(&data, ref u1b, ref u2b,
                out _, out _, out double rhoB, out _,
                1e-14, out _, 30, g2, H2);

            // Both should converge to critical points with the same rho
            AssertClose(rhoA, rhoB, 1e-6, "symmetric rho");
        }
    }
}
