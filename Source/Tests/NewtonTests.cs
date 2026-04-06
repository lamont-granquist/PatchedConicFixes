using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public unsafe class NewtonTests
    {
        private const double Tol = 1e-12;

        private static void BuildPolynomial(Cmplx[] roots, Cmplx* c)
        {
            int n = roots.Length;

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

        private static Cmplx EvalPoly(Cmplx* c, int n, Cmplx z)
        {
            Cmplx val = c[n];
            for (int i = n - 1; i >= 0; i--)
                val = val * z + c[i];
            return val;
        }

        private static void AssertIsRoot(Cmplx* c, int n, Cmplx root, double tol = Tol)
        {
            Cmplx  val   = EvalPoly(c, n, root);
            double scale = Math.Max(1.0, Math.Pow(root.Abs, n));
            Assert.True(val.Abs < tol * scale,
                $"|P({root.Re:G17} + {root.Im:G17}i)| = {val.Abs:G17} exceeds {tol * scale:G17}");
        }

        private static void AssertNear(Cmplx expected, Cmplx actual, double tol = Tol)
        {
            Assert.True((expected - actual).Abs < tol,
                $"Expected ({expected.Re:G17}, {expected.Im:G17}) " +
                $"got ({actual.Re:G17}, {actual.Im:G17}), " +
                $"diff = {(expected - actual).Abs:G17}");
        }

        // --- basic convergence to known roots ---

        [Fact]
        public void SimpleRealRoot_CloseStart()
        {
            // (z - 3)(z + 1) = z² - 2z - 3
            Cmplx* c     = stackalloc Cmplx[3];
            var    roots = new Cmplx[] { 3.0, -1.0 };
            BuildPolynomial(roots, c);

            Cmplx z = 2.8;
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(3.0, z);
            AssertIsRoot(c, 2, z);
        }

        [Fact]
        public void SimpleRealRoot_FarStart()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            var    roots = new Cmplx[] { 3.0, -1.0 };
            BuildPolynomial(roots, c);

            Cmplx z = 100.0;
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertIsRoot(c, 2, z);
        }

        [Fact]
        public void NegativeRealRoot()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            var    roots = new Cmplx[] { 3.0, -1.0 };
            BuildPolynomial(roots, c);

            Cmplx z = -0.8;
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(-1.0, z);
            AssertIsRoot(c, 2, z);
        }

        // --- complex roots ---

        [Fact]
        public void ComplexConjugateRoot_FromNearby()
        {
            // z² + 1 = 0, roots ±i
            Cmplx*  c     = stackalloc Cmplx[3];
            Cmplx[] roots = { new Cmplx(0, 1), new Cmplx(0, -1) };
            BuildPolynomial(roots, c);

            var z = new Cmplx(0.1, 1.1);
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(new Cmplx(0, 1), z);
        }

        [Fact]
        public void ComplexRoot_GeneralPosition()
        {
            Cmplx*  c     = stackalloc Cmplx[4];
            var     r0    = new Cmplx(1.0, 2.0);
            var     r1    = new Cmplx(1.0, -2.0);
            var     r2    = new Cmplx(-3.0);
            Cmplx[] roots = { r0, r1, r2 };
            BuildPolynomial(roots, c);

            var z = new Cmplx(0.8, 2.2);
            PolySolver.Newton(3, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(r0, z);
            AssertIsRoot(c, 3, z);
        }

        // --- roots on the unit circle (MOID-like) ---

        [Fact]
        public void RootOnUnitCircle_FromNearby()
        {
            Cmplx* c  = stackalloc Cmplx[3];
            var    z0 = Cmplx.Polar(1.0, Math.PI / 3.0);
            var    z1 = Cmplx.Polar(1.0, -Math.PI / 3.0);
            BuildPolynomial(new[] { z0, z1 }, c);

            var z = Cmplx.Polar(1.05, Math.PI / 3.0 + 0.05);
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(z0, z);
        }

        // --- root inside unit circle (inner path) ---

        [Fact]
        public void RootInsideUnitCircle()
        {
            Cmplx* c  = stackalloc Cmplx[3];
            var    r0 = new Cmplx(0.3, 0.4); // |z| = 0.5
            var    r1 = new Cmplx(-0.1, -0.2);
            BuildPolynomial(new[] { r0, r1 }, c);

            var z = new Cmplx(0.35, 0.45);
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(r0, z);
        }

        // --- root well outside unit circle (outer path) ---

        [Fact]
        public void RootOutsideUnitCircle()
        {
            Cmplx* c  = stackalloc Cmplx[3];
            var    r0 = new Cmplx(5.0, 3.0); // |z| ≈ 5.83
            var    r1 = new Cmplx(0.2);
            BuildPolynomial(new[] { r0, r1 }, c);

            var z = new Cmplx(5.5, 3.5);
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(r0, z);
        }

        // --- root at the inner/outer switching boundary ---

        [Fact]
        public void RootNearSwitchingThreshold()
        {
            // |z|² ≈ 5 is the switching boundary
            Cmplx* c  = stackalloc Cmplx[3];
            var    r0 = Cmplx.Polar(Math.Sqrt(5.0), 0.7); // exactly at boundary
            var    r1 = new Cmplx(0.1);
            BuildPolynomial(new[] { r0, r1 }, c);

            Cmplx z = r0 + new Cmplx(0.01, -0.01);
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(r0, z);
        }

        // --- near-double root (ill-conditioned) ---

        [Fact]
        public void NearDoubleRoot()
        {
            Cmplx* c  = stackalloc Cmplx[3];
            var    r0 = new Cmplx(2.0);
            var    r1 = new Cmplx(2.0 + 1e-6);
            BuildPolynomial(new[] { r0, r1 }, c);

            Cmplx z = 1.9;
            PolySolver.Newton(2, c, ref z, 1e-8, 0.0, 300, 10, PolySolver.MakeRng());

            // Relaxed tolerance — near-double roots are ill-conditioned
            AssertIsRoot(c, 2, z, 1e-8);
        }

        [Fact]
        public void ExactDoubleRoot()
        {
            // (z - 4)²(z + 1) — Newton should still converge, just slowly
            Cmplx* c = stackalloc Cmplx[4];
            BuildPolynomial(new Cmplx[] { 4.0, 4.0, -1.0 }, c);

            Cmplx z = 3.5;
            PolySolver.Newton(3, c, ref z, 1e-6, 0.0, 300, 10, PolySolver.MakeRng());

            AssertIsRoot(c, 3, z, 1e-6);
        }

        // --- starting from zero ---

        [Fact]
        public void StartFromZero()
        {
            Cmplx* c = stackalloc Cmplx[3];
            BuildPolynomial(new Cmplx[] { 2.0, -3.0 }, c);

            Cmplx z = 0.0;
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertIsRoot(c, 2, z);
        }

        // --- starting from a root directly ---

        [Fact]
        public void StartAtExactRoot()
        {
            Cmplx* c = stackalloc Cmplx[3];
            BuildPolynomial(new Cmplx[] { 2.0, -3.0 }, c);

            Cmplx z     = 2.0;
            ulong iters = PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(2.0, z);
            // Should converge immediately
            Assert.True(iters <= 3, $"Expected few iterations at exact root, got {iters}");
        }

        // --- iteration count is reasonable ---

        [Fact]
        public void IterationCountIsReasonable_SimpleCase()
        {
            Cmplx* c = stackalloc Cmplx[3];
            BuildPolynomial(new Cmplx[] { 1.0, -1.0 }, c);

            Cmplx z     = 1.5;
            ulong iters = PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(1.0, z);
            // Newton on a quadratic from nearby should converge in well under 50
            Assert.True(iters < 50, $"Too many iterations: {iters}");
        }

        [Fact]
        public void IterationCountIsReasonable_Degree16()
        {
            const int N     = 16;
            Cmplx*    c     = stackalloc Cmplx[N + 1];
            var       roots = new Cmplx[N];

            // Dense polynomial with roots near |z| = 1, mimics MOID structure
            for (int i = 0; i < N; i++)
                roots[i] = Cmplx.Polar(1.0 + 0.01 * Math.Sin(i), 2.0 * Math.PI * i / N);
            BuildPolynomial(roots, c);

            // Start near roots[1]
            Cmplx z     = roots[1] * Cmplx.Polar(1.0, 0.02);
            ulong iters = PolySolver.Newton(N, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertIsRoot(c, N, z);
            Assert.True(iters < 100, $"Too many iterations for degree 16: {iters}");
        }

        [Fact]
        public void RootOnUnitCircle_Degree8()
        {
            const int N     = 8;
            Cmplx*    c     = stackalloc Cmplx[N + 1];
            var       roots = new Cmplx[N];

            // Dense polynomial with roots on the unit circle
            for (int i = 0; i < N; i++)
                roots[i] = Cmplx.Polar(1.0, 2.0 * Math.PI * i / N + 0.03 * i);
            BuildPolynomial(roots, c);

            Cmplx target = roots[1];
            Cmplx z      = target * Cmplx.Polar(1.0, 0.02);
            PolySolver.Newton(N, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertIsRoot(c, N, z);
        }

        // --- cycle detection exercises ---

        [Fact]
        public void ConvergesFromPoorStartingGuess()
        {
            // Degree 5 polynomial, starting far from all roots
            Cmplx*  c     = stackalloc Cmplx[6];
            Cmplx[] roots = { new Cmplx(1, 1), new Cmplx(1, -1), new Cmplx(-1, 1), new Cmplx(-1, -1), new Cmplx(0.5) };
            BuildPolynomial(roots, c);

            var z = new Cmplx(50.0, 50.0);
            PolySolver.Newton(5, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertIsRoot(c, 5, z);
        }

        [Fact]
        public void ConvergesFromNegativeRealAxis()
        {
            // Start on the negative real axis, far from complex roots
            Cmplx*  c     = stackalloc Cmplx[5];
            Cmplx[] roots = { Cmplx.Polar(1.0, 0.5), Cmplx.Polar(1.0, -0.5), Cmplx.Polar(1.0, 2.0), Cmplx.Polar(1.0, -2.0) };
            BuildPolynomial(roots, c);

            Cmplx z = -10.0;
            PolySolver.Newton(4, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertIsRoot(c, 4, z);
        }

        // --- mineps = 0 (squeeze to machine precision) ---

        [Fact]
        public void MinEpsZero_ConvergesToMachinePrecision()
        {
            Cmplx* c = stackalloc Cmplx[3];
            BuildPolynomial(new Cmplx[] { 2.0, -3.0 }, c);

            Cmplx z = 1.5;
            PolySolver.Newton(2, c, ref z, 1e-8, 0.0, 300, 10, PolySolver.MakeRng());

            // Should be extremely close to 2.0
            AssertNear(2.0, z, 1e-14);
        }

        // --- root at zero ---

        [Fact]
        public void RootAtZero()
        {
            // z(z - 2) = z² - 2z
            Cmplx* c = stackalloc Cmplx[3];
            BuildPolynomial(new Cmplx[] { 0.0, 2.0 }, c);

            Cmplx z = 0.01;
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertIsRoot(c, 2, z);
        }

        // --- MOID-like degree-16 trig polynomial structure ---

        // NOTE: This solver isn't robust against sparse cases like z^16 - 1 = 0, but the MOID polynomial
        // structure is never that ill-conditioned.
        [Fact]
        public void Degree16_TrigPolynomialRoots()
        {
            // Build a polynomial with conjugate-symmetric structure like the MOID polynomial:
            // roots come in pairs (z, 1/z*), all near |z| = 1
            const int N     = 16;
            Cmplx*    c     = stackalloc Cmplx[N + 1];
            var       roots = new Cmplx[N];

            // 8 pairs with slight perturbation off unit circle
            for (int i = 0; i < 8; i++)
            {
                double angle = Math.PI * i / 8.0 + 0.1;
                double r     = 1.0 + 0.01 * (i % 3 - 1); // 0.99, 1.0, or 1.01
                roots[2 * i]     = Cmplx.Polar(r, angle);
                roots[2 * i + 1] = Cmplx.Polar(1.0 / r, -angle); // 1/z*
            }

            BuildPolynomial(roots, c);

            // Try converging to a few of them from nearby starts
            for (int i = 0; i < 4; i++)
            {
                Cmplx z = roots[i] * Cmplx.Polar(1.0, 0.05); // slight angular perturbation
                PolySolver.Newton(N, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());
                AssertIsRoot(c, N, z);
            }
        }

        // --- non-monic polynomial ---

        [Fact]
        public void NonMonicPolynomial()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            var    roots = new Cmplx[] { 1.0, -2.0 };
            BuildPolynomial(roots, c);
            var scale                        = new Cmplx(7.0, -3.0);
            for (int i = 0; i < 3; i++) c[i] = c[i] * scale;

            Cmplx z = 0.8;
            PolySolver.Newton(2, c, ref z, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            AssertNear(1.0, z);
            AssertIsRoot(c, 2, z);
        }

        // --- multiple calls find different roots with different starts ---

        [Fact]
        public void DifferentStartsCanFindDifferentRoots()
        {
            Cmplx*  c     = stackalloc Cmplx[4];
            Cmplx[] roots = { new Cmplx(1.0), new Cmplx(-2.0), new Cmplx(0.0, 3.0) };
            BuildPolynomial(roots, c);

            Cmplx z1 = 0.9;
            PolySolver.Newton(3, c, ref z1, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            Cmplx z2 = -1.8;
            PolySolver.Newton(3, c, ref z2, 1e-8, 1e-14, 300, 10, PolySolver.MakeRng());

            // Both should be valid roots
            AssertIsRoot(c, 3, z1);
            AssertIsRoot(c, 3, z2);

            // But they should be different roots
            Assert.True((z1 - z2).Abs > 0.5,
                $"Expected different roots but got z1=({z1.Re:G4},{z1.Im:G4}) z2=({z2.Re:G4},{z2.Im:G4})");
        }
    }
}
