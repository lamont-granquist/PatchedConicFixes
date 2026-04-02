using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public unsafe class RootBoundAndErrorTests
    {
        private const double Tol = 1e-10;

        private static void AssertClose(double expected, double actual, double tol = Tol, string label = "")
        {
            Assert.True(Math.Abs(expected - actual) < tol,
                $"{label} Expected {expected:G17}, got {actual:G17}, diff = {Math.Abs(expected - actual):G17}");
        }

        private static void BuildPolynomial(Cmplx[] roots, Cmplx* c)
        {
            int n                             = roots.Length;
            for (int i = 0; i <= n; i++) c[i] = 0.0;
            c[0] = 1.0;
            for (int k = 0; k < n; k++)
            {
                c[k + 1] = c[k];
                for (int i = k; i > 0; i--)
                    c[i] = c[i - 1] - roots[k] * c[i];
                c[0] = -roots[k] * c[0];
            }
        }

        // =============================================================
        // MaxRootBound tests
        // =============================================================

        [Fact]
        public void MaxRootBound_Forward_ReasonableEstimate()
        {
            Cmplx*  c     = stackalloc Cmplx[5];
            Cmplx[] roots = { new Cmplx(0.5, 0.3), new Cmplx(-1.2, 0.8), new Cmplx(2.0), new Cmplx(-0.1, -0.9) };
            BuildPolynomial(roots, c);

            double bound = Baluev.MaxRootBound(4, c, true);

            double maxActual = 0;
            foreach (Cmplx z in roots)
                if (z.Abs > maxActual)
                    maxActual = z.Abs;

            // Should be in the right ballpark — within a small factor of the actual max
            Assert.True(bound > 0, "Bound should be positive");
            Assert.True(bound < maxActual * 10,
                $"Bound {bound:G17} is unreasonably large vs max |root| = {maxActual:G17}");
        }

        [Fact]
        public void MaxRootBound_Backward_ReasonableEstimate()
        {
            Cmplx*  c     = stackalloc Cmplx[5];
            Cmplx[] roots = { new Cmplx(0.5, 0.3), new Cmplx(-1.2, 0.8), new Cmplx(2.0), new Cmplx(-0.1, -0.9) };
            BuildPolynomial(roots, c);

            double bound = Baluev.MaxRootBound(4, c, false);

            Assert.True(bound > 0, "Bound should be positive");
            Assert.True(bound < 100, "Backward bound should be reasonable");
        }

        [Fact]
        public void MaxRootBound_UnitCircleRoots()
        {
            const int N     = 8;
            Cmplx*    c     = stackalloc Cmplx[N + 1];
            var       roots = new Cmplx[N];
            for (int i = 0; i < N; i++)
                roots[i] = Cmplx.Polar(1.0, 2.0 * Math.PI * i / N + 0.05 * i);
            BuildPolynomial(roots, c);

            double bound = Baluev.MaxRootBound(N, c, true);

            // For unit circle roots the estimate should be close to 1
            Assert.True(bound > 0.5 && bound < 5.0,
                $"Bound {bound:G17} should be near 1.0 for unit circle roots");
        }

        [Fact]
        public void MaxRootBound_Degree16_MOIDLike()
        {
            const int N     = 16;
            Cmplx*    c     = stackalloc Cmplx[N + 1];
            var       roots = new Cmplx[N];
            for (int i = 0; i < 8; i++)
            {
                double angle = Math.PI * i / 8.0 + 0.15;
                double r     = 1.0 + 0.02 * Math.Sin(i);
                roots[2 * i]     = Cmplx.Polar(r, angle);
                roots[2 * i + 1] = Cmplx.Polar(1.0 / r, -angle);
            }

            BuildPolynomial(roots, c);

            double forward  = Baluev.MaxRootBound(N, c, true);
            double backward = Baluev.MaxRootBound(N, c, false);

            // For MOID-like roots near |z|=1, both bounds should be close to 1
            Assert.True(forward > 0.5 && forward < 5.0,
                $"Forward bound {forward:G17} should be near 1.0");
            Assert.True(backward > 0.5 && backward < 5.0,
                $"Backward bound {backward:G17} should be near 1.0");
        }

        [Fact]
        public void MaxRootBound_ForwardGteBackward_ForMixedRoots()
        {
            // Forward bound estimates max |z|, backward estimates max |1/z|
            // For roots both inside and outside unit circle, forward >= backward isn't guaranteed
            // but we can check they're both positive and finite
            Cmplx*  c     = stackalloc Cmplx[5];
            Cmplx[] roots = { new Cmplx(0.5, 0.3), new Cmplx(-1.2, 0.8), new Cmplx(2.0), new Cmplx(-0.1, -0.9) };
            BuildPolynomial(roots, c);

            double fwd = Baluev.MaxRootBound(4, c, true);
            double bwd = Baluev.MaxRootBound(4, c, false);

            Assert.True(fwd > 0 && !double.IsInfinity(fwd));
            Assert.True(bwd > 0 && !double.IsInfinity(bwd));
        }

        [Fact]
        public void MaxRootBound_SmallDegree()
        {
            // Degree 1: z - 3 = 0, root = 3
            Cmplx* c = stackalloc Cmplx[2];
            c[0] = -3.0;
            c[1] = 1.0;

            double bound = Baluev.MaxRootBound(1, c, true);
            Assert.True(bound > 0, $"Bound should be positive, got {bound:G17}");
        }

        [Fact]
        public void MaxRootBound_DegreeZero_ReturnsInfinity()
        {
            Cmplx* c = stackalloc Cmplx[1];
            c[0] = 5.0;

            double bound = Baluev.MaxRootBound(0, c, true);
            Assert.True(double.IsPositiveInfinity(bound));
        }

        [Fact]
        public void MaxRootBound_NonMonic()
        {
            Cmplx* c     = stackalloc Cmplx[4];
            var    roots = new Cmplx[] { 1.0, -2.0, 3.0 };
            BuildPolynomial(roots, c);
            var scale                         = new Cmplx(5.0, -2.0);
            for (int i = 0; i <= 3; i++) c[i] = c[i] * scale;

            double bound = Baluev.MaxRootBound(3, c, true);

            // Should be positive, finite, and in the right ballpark
            Assert.True(bound > 0 && !double.IsInfinity(bound),
                $"Bound should be positive and finite, got {bound:G17}");
        }


        // =============================================================
        // RootError tests
        // =============================================================

        [Fact]
        public void RootError_ExactRoot_SmallError()
        {
            // Build polynomial from known roots, evaluate error at an exact root
            Cmplx*  c     = stackalloc Cmplx[5];
            Cmplx[] roots = { new Cmplx(1.0, 2.0), new Cmplx(1.0, -2.0), new Cmplx(-1.0, 0.5), new Cmplx(-1.0, -0.5) };
            BuildPolynomial(roots, c);

            double err = Baluev.RootError(4, c, 1e-15, roots[0]);

            Assert.True(err < 1e-10,
                $"Error at exact root should be tiny, got {err:G17}");
        }

        [Fact]
        public void RootError_ScalesWithCerr()
        {
            // Larger cerr should produce larger root error
            Cmplx* c     = stackalloc Cmplx[4];
            var    roots = new Cmplx[] { 1.0, -2.0, 3.0 };
            BuildPolynomial(roots, c);

            double errSmall = Baluev.RootError(3, c, 1e-15, roots[0]);
            double errLarge = Baluev.RootError(3, c, 1e-8, roots[0]);

            Assert.True(errLarge > errSmall,
                $"Larger cerr should give larger error: {errSmall:G17} vs {errLarge:G17}");
        }

        [Fact]
        public void RootError_IsPositive()
        {
            Cmplx*  c     = stackalloc Cmplx[4];
            Cmplx[] roots = { new Cmplx(0.5, 0.5), new Cmplx(0.5, -0.5), new Cmplx(2.0) };
            BuildPolynomial(roots, c);

            foreach (Cmplx z in roots)
            {
                double err = Baluev.RootError(3, c, 1e-14, z);
                Assert.True(err >= 0, $"Error should be non-negative at {z}, got {err:G17}");
            }
        }

        [Fact]
        public void RootError_InsideUnitCircle_UsesForwardPath()
        {
            // Root with |z| < 1 uses forward evaluation
            Cmplx* c = stackalloc Cmplx[3];
            Cmplx[] roots =
            {
                new Cmplx(0.3, 0.4), // |z| = 0.5
                new Cmplx(0.3, -0.4)
            };
            BuildPolynomial(roots, c);

            double err = Baluev.RootError(2, c, 1e-14, roots[0]);
            Assert.True(err < 1e-8 && err >= 0,
                $"Error at exact root inside unit circle: {err:G17}");
        }

        [Fact]
        public void RootError_OutsideUnitCircle_UsesBackwardPath()
        {
            // Root with |z| > 1 uses backward evaluation
            Cmplx* c = stackalloc Cmplx[3];
            Cmplx[] roots =
            {
                new Cmplx(3.0, 4.0), // |z| = 5
                new Cmplx(3.0, -4.0)
            };
            BuildPolynomial(roots, c);

            double err = Baluev.RootError(2, c, 1e-14, roots[0]);
            Assert.True(err < 1e-8 && err >= 0,
                $"Error at exact root outside unit circle: {err:G17}");
        }

        [Fact]
        public void RootError_OnUnitCircle()
        {
            // Root with |z| = 1 (boundary of forward/backward)
            Cmplx*  c     = stackalloc Cmplx[3];
            var     z0    = Cmplx.Polar(1.0, 0.7);
            var     z1    = Cmplx.Polar(1.0, -0.7);
            Cmplx[] roots = { z0, z1 };
            BuildPolynomial(roots, c);

            double err = Baluev.RootError(2, c, 1e-14, z0);
            Assert.True(err < 1e-8 && err >= 0,
                $"Error at unit circle root: {err:G17}");
        }

        [Fact]
        public void RootError_NearDoubleRoot_LargerError()
        {
            // Near-double roots are ill-conditioned, error should be larger
            Cmplx*  cWell   = stackalloc Cmplx[3];
            Cmplx[] wellSep = { new Cmplx(1.0), new Cmplx(-1.0) };
            BuildPolynomial(wellSep, cWell);

            Cmplx*  cNear      = stackalloc Cmplx[3];
            Cmplx[] nearDouble = { new Cmplx(1.0), new Cmplx(1.0 + 1e-6) };
            BuildPolynomial(nearDouble, cNear);

            double errWell = Baluev.RootError(2, cWell, 1e-14, wellSep[0]);
            double errNear = Baluev.RootError(2, cNear, 1e-14, nearDouble[0]);

            Assert.True(errNear > errWell,
                $"Near-double root should have larger error: well={errWell:G17} near={errNear:G17}");
        }

        [Fact]
        public void RootError_Degree16_MOIDLike()
        {
            const int N     = 16;
            Cmplx*    c     = stackalloc Cmplx[N + 1];
            var       roots = new Cmplx[N];
            for (int i = 0; i < 8; i++)
            {
                double angle = Math.PI * i / 8.0 + 0.15;
                double r     = 1.0 + 0.02 * Math.Sin(i);
                roots[2 * i]     = Cmplx.Polar(r, angle);
                roots[2 * i + 1] = Cmplx.Polar(1.0 / r, -angle);
            }

            BuildPolynomial(roots, c);

            // All roots should have small error with small cerr
            foreach (Cmplx z in roots)
            {
                double err = Baluev.RootError(N, c, 1e-14, z);
                Assert.True(err < 1e-6,
                    $"Error at MOID-like root {z} = {err:G17}");
            }
        }

        [Fact]
        public void RootError_PerturbedRoot_ErrorBoundsActualError()
        {
            // Perturb a root slightly and check that RootError gives
            // a reasonable bound on the actual perturbation
            Cmplx* c     = stackalloc Cmplx[4];
            var    roots = new Cmplx[] { 1.0, -2.0, 3.0 };
            BuildPolynomial(roots, c);

            double perturbation = 1e-8;
            Cmplx  perturbed    = roots[0] + perturbation;

            double estimatedErr = Baluev.RootError(3, c, 1e-15, perturbed);

            // The estimated error should be in the right ballpark of the actual perturbation
            // (not orders of magnitude off in either direction)
            Assert.True(estimatedErr > perturbation * 0.01,
                $"Estimated error {estimatedErr:G17} is unreasonably small for perturbation {perturbation}");
            Assert.True(estimatedErr < perturbation * 1000,
                $"Estimated error {estimatedErr:G17} is unreasonably large for perturbation {perturbation}");
        }

        [Fact]
        public void RootError_MixedMagnitudes()
        {
            // Roots both inside and outside unit circle
            Cmplx* c = stackalloc Cmplx[5];
            Cmplx[] roots =
            {
                new Cmplx(0.2, 0.1),  // inside
                new Cmplx(0.2, -0.1), // inside
                new Cmplx(3.0, 1.0),  // outside
                new Cmplx(3.0, -1.0)  // outside
            };
            BuildPolynomial(roots, c);

            foreach (Cmplx z in roots)
            {
                double err = Baluev.RootError(4, c, 1e-14, z);
                Assert.True(err >= 0 && err < 1e-6,
                    $"Error at root {z} = {err:G17}");
            }
        }
    }
}
