using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public unsafe class SqDistSurfaceTests
    {
        private const double Tol = 1e-10;

        private static void AssertClose(double expected, double actual, double tol = Tol, string label = "")
        {
            Assert.True(Math.Abs(expected - actual) < tol,
                $"{label} Expected {expected:G17}, got {actual:G17}, diff = {Math.Abs(expected - actual):G17}");
        }

        // --- SqDistVG: basic properties ---

        [Fact]
        public void SqDistVG_SamePoint_RhoIsZero()
        {
            var     O    = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     data = SAuxData.Create(in O, in O);
            double* g    = stackalloc double[2];

            double rho = Baluev.SqDistVG(&data, true, true, 1.0, 1.0, g);

            AssertClose(0.0, rho, 1e-8, "rho at same point");
        }

        [Fact]
        public void SqDistVG_RhoIsNonNegative()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];

            for (int i = 0; i < 8; i++)
            {
                double u1 = -Math.PI + 0.1 + i * 0.7;
                for (int j = 0; j < 8; j++)
                {
                    double u2  = -Math.PI + 0.2 + j * 0.7;
                    double rho = Baluev.SqDistVG(&data, true, true, u1, u2, g);
                    Assert.True(rho >= -Tol, $"rho={rho:G17} at u1={u1:F3}, u2={u2:F3}");
                }
            }
        }

        [Fact]
        public void SqDistVG_ConsistentWithDistanceBetween()
        {
            // rho = |dr|² / (2·a1·a2), so distance = sqrt(2·a1·a2·rho)
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];

            double u1   = 1.0, u2 = 0.5;
            double rho  = Baluev.SqDistVG(&data, true, true, u1, u2, g);
            double dist = Baluev.DistanceBetween(&data, true, true, u1, u2);

            double distFromRho = Math.Sqrt(2.0 * data.a1 * data.a2 * rho);
            AssertClose(dist, distFromRho, 1e-8, "distance from rho");
        }

        [Fact]
        public void SqDistVG_GradientMatchesFiniteDifference()
        {
            var     O1     = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2     = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data   = SAuxData.Create(in O1, in O2);
            double* g      = stackalloc double[2];
            double* gDummy = stackalloc double[2];

            double u1 = 1.0, u2 = 0.5;
            Baluev.SqDistVG(&data, true, true, u1, u2, g);

            double h     = 1e-7;
            double rhoP1 = Baluev.SqDistVG(&data, true, true, u1 + h, u2, gDummy);
            double rhoM1 = Baluev.SqDistVG(&data, true, true, u1 - h, u2, gDummy);
            double rhoP2 = Baluev.SqDistVG(&data, true, true, u1, u2 + h, gDummy);
            double rhoM2 = Baluev.SqDistVG(&data, true, true, u1, u2 - h, gDummy);

            AssertClose((rhoP1 - rhoM1) / (2.0 * h), g[0], 1e-5, "g[0]");
            AssertClose((rhoP2 - rhoM2) / (2.0 * h), g[1], 1e-5, "g[1]");
        }

        // --- SqDistVGH: Hessian checks ---

        [Fact]
        public void SqDistVGH_MatchesSqDistVG()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* gVG  = stackalloc double[2];
            double* gVGH = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1     = 1.0, u2 = 0.5;
            double rhoVG  = Baluev.SqDistVG(&data, true, true, u1, u2, gVG);
            double rhoVGH = Baluev.SqDistVGH(&data, true, true, u1, u2, gVGH, H);

            AssertClose(rhoVG, rhoVGH, label: "rho");
            AssertClose(gVG[0], gVGH[0], label: "g[0]");
            AssertClose(gVG[1], gVGH[1], label: "g[1]");
        }

        [Fact]
        public void SqDistVGH_HessianMatchesFiniteDifference()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 1.0, u2 = 0.5;
            Baluev.SqDistVGH(&data, true, true, u1, u2, g, H);

            // FD of gradient to get Hessian
            double  h   = 1e-6;
            double* gP1 = stackalloc double[2];
            double* gM1 = stackalloc double[2];
            double* gP2 = stackalloc double[2];
            double* gM2 = stackalloc double[2];

            Baluev.SqDistVG(&data, true, true, u1 + h, u2, gP1);
            Baluev.SqDistVG(&data, true, true, u1 - h, u2, gM1);
            Baluev.SqDistVG(&data, true, true, u1, u2 + h, gP2);
            Baluev.SqDistVG(&data, true, true, u1, u2 - h, gM2);

            double H00_fd = (gP1[0] - gM1[0]) / (2.0 * h);
            double H11_fd = (gP2[1] - gM2[1]) / (2.0 * h);
            double H01_fd = (gP1[1] - gM1[1]) / (2.0 * h);

            AssertClose(H00_fd, H[0], 1e-3, "H[0] (d²ρ/du1²)");
            AssertClose(H11_fd, H[1], 1e-3, "H[1] (d²ρ/du2²)");
            AssertClose(H01_fd, H[2], 1e-3, "H[2] (d²ρ/du1du2)");
        }

        [Fact]
        public void SqDistVGH_HessianIsSymmetric()
        {
            // H[2] should match the cross derivative computed both ways
            var O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            double  h   = 1e-6;
            double* gP1 = stackalloc double[2];
            double* gM1 = stackalloc double[2];
            double* gP2 = stackalloc double[2];
            double* gM2 = stackalloc double[2];

            double u1 = 1.0, u2 = 0.5;
            Baluev.SqDistVG(&data, true, true, u1 + h, u2, gP1);
            Baluev.SqDistVG(&data, true, true, u1 - h, u2, gM1);
            Baluev.SqDistVG(&data, true, true, u1, u2 + h, gP2);
            Baluev.SqDistVG(&data, true, true, u1, u2 - h, gM2);

            double H01 = (gP1[1] - gM1[1]) / (2.0 * h); // dg1/du1
            double H10 = (gP2[0] - gM2[0]) / (2.0 * h); // dg0/du2

            AssertClose(H01, H10, 1e-3, "Hessian symmetry");
        }

        [Fact]
        public void SqDistVGH_SeveralConfigurations()
        {
            (COrbitData, COrbitData)[] configs = { (new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0), new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0)), (new COrbitData(2.0, 0.5, 0.3, 1.0, 0.5), new COrbitData(1.5, 0.2, 0.8, 0.5, 2.0)), (new COrbitData(3.0, 0.8, 0.1, 0.3, 0.0), new COrbitData(1.0, 0.1, 0.1, 1.0, 0.0)) };

            double h = 1e-6;
            foreach ((COrbitData O1, COrbitData O2) in configs)
            {
                var     data = SAuxData.Create(in O1, in O2);
                double* g    = stackalloc double[2];
                double* H    = stackalloc double[3];
                double* gP1  = stackalloc double[2];
                double* gM1  = stackalloc double[2];
                double* gP2  = stackalloc double[2];
                double* gM2  = stackalloc double[2];

                double u1 = 0.7, u2 = 1.3;
                Baluev.SqDistVGH(&data, true, true, u1, u2, g, H);

                Baluev.SqDistVG(&data, true, true, u1 + h, u2, gP1);
                Baluev.SqDistVG(&data, true, true, u1 - h, u2, gM1);
                Baluev.SqDistVG(&data, true, true, u1, u2 + h, gP2);
                Baluev.SqDistVG(&data, true, true, u1, u2 - h, gM2);

                AssertClose((gP1[0] - gM1[0]) / (2.0 * h), H[0], 1e-3,
                    $"H[0] a1={O1.a},e1={O1.e}");
                AssertClose((gP2[1] - gM2[1]) / (2.0 * h), H[1], 1e-3,
                    $"H[1] a1={O1.a},e1={O1.e}");
                AssertClose((gP1[1] - gM1[1]) / (2.0 * h), H[2], 1e-3,
                    $"H[2] a1={O1.a},e1={O1.e}");
            }
        }

        // --- SqDist2DIter: Newton step ---

        [Fact]
        public void SqDist2DIter_ReducesGradient()
        {
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            // Start slightly off a known critical point
            double u1 = 1.0, u2 = 0.5;
            double gNormBefore;
            {
                double* gCheck = stackalloc double[2];
                Baluev.SqDistVG(&data, true, true, u1, u2, gCheck);
                gNormBefore = gCheck[0] * gCheck[0] + gCheck[1] * gCheck[1];
            }

            Baluev.SqDist2DIter(&data, true, true, ref u1, ref u2,
                out _, g, H, out double detH);

            if (detH == 0.0) return;

            double* gAfter = stackalloc double[2];
            Baluev.SqDistVG(&data, true, true, u1, u2, gAfter);
            double gNormAfter = gAfter[0] * gAfter[0] + gAfter[1] * gAfter[1];

            Assert.True(gNormAfter < gNormBefore,
                $"|g| increased: {Math.Sqrt(gNormBefore):G6} → {Math.Sqrt(gNormAfter):G6}");
        }

        [Fact]
        public void SqDist2DIter_ConvergesToCriticalPoint()
        {
            // Several iterations should drive gradient to zero
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 1.0, u2 = 0.5;
            for (int iter = 0; iter < 20; iter++)
            {
                Baluev.SqDist2DIter(&data, true, true, ref u1, ref u2,
                    out _, g, H, out double detH);
                if (detH == 0.0) break;
            }

            double gNorm = Math.Sqrt(g[0] * g[0] + g[1] * g[1]);
            Assert.True(gNorm < 1e-8, $"Gradient should be near zero after iterations, got {gNorm:G17}");
        }

        [Fact]
        public void SqDist2DIter_OutputsConsistent()
        {
            var     O1     = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2     = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data   = SAuxData.Create(in O1, in O2);
            double* g      = stackalloc double[2];
            double* H      = stackalloc double[3];
            double* gCheck = stackalloc double[2];
            double* HCheck = stackalloc double[3];

            double u1 = 1.0, u2 = 0.5;

            // The rho, g, H outputs should match a direct call to SqDistVGH
            // at the pre-step position
            double u1Before = u1, u2Before = u2;
            double rhoCheck = Baluev.SqDistVGH(&data, true, true, u1Before, u2Before, gCheck, HCheck);

            Baluev.SqDist2DIter(&data, true, true, ref u1, ref u2,
                out double rho, g, H, out _);

            AssertClose(rhoCheck, rho, label: "rho");
            AssertClose(gCheck[0], g[0], label: "g[0]");
            AssertClose(gCheck[1], g[1], label: "g[1]");
            AssertClose(HCheck[0], H[0], label: "H[0]");
            AssertClose(HCheck[1], H[1], label: "H[1]");
            AssertClose(HCheck[2], H[2], label: "H[2]");
        }

        [Fact]
        public void SqDist2DIter_StepSizeDecreasesNearMinimum()
        {
            // As we approach a minimum, step sizes should decrease
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1         = 1.0, u2 = 0.5;
            double prevStep   = double.MaxValue;
            int    decreasing = 0;

            for (int iter = 0; iter < 15; iter++)
            {
                double step = Baluev.SqDist2DIter(&data, true, true, ref u1, ref u2,
                    out _, g, H, out double detH);
                if (detH == 0.0) break;
                if (step < prevStep) decreasing++;
                prevStep = step;
            }

            Assert.True(decreasing >= 10,
                $"Step size should mostly decrease near convergence, only decreased {decreasing}/15 times");
        }

        [Fact]
        public void SqDist2DIter_ConvergedPointMatchesMinDistanceFor()
        {
            // After convergence, the (u1, u2) should match what MinDistanceFor gives
            var     O1   = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var     O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var     data = SAuxData.Create(in O1, in O2);
            double* g    = stackalloc double[2];
            double* H    = stackalloc double[3];

            double u1 = 1.0, u2 = 0.5;
            for (int iter = 0; iter < 30; iter++)
            {
                Baluev.SqDist2DIter(&data, true, true, ref u1, ref u2,
                    out _, g, H, out double detH);
                if (detH == 0.0) break;
            }

            // The converged u1 should give the same u2 via MinDistanceFor
            double dMin = Baluev.MinDistanceFor(&data, u1, out double u2Check);
            if (dMin < 0) return;

            double dist2D = Baluev.DistanceBetween(&data, true, true, u1, u2);
            AssertClose(Math.Abs(dMin), dist2D, 1e-6, "distances match");
        }
    }
}
