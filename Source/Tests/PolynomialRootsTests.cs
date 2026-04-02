using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public unsafe class PolynomialRootsTests
    {
        private const double Tol = 1e-10;

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

        private static Cmplx EvalPoly(Cmplx* c, int n, Cmplx z)
        {
            Cmplx val = c[n];
            for (int i = n - 1; i >= 0; i--)
                val = val * z + c[i];
            return val;
        }

        private static void CopyCoeffs(Cmplx* src, Cmplx* dst, int n)
        {
            for (int i = 0; i <= n; i++) dst[i] = src[i];
        }

        private static void AssertAllResiduals(Cmplx* c, int n, Cmplx* roots, double tol = Tol)
        {
            for (int i = 0; i < n; i++)
            {
                Cmplx  val   = EvalPoly(c, n, roots[i]);
                double scale = Math.Max(1.0, Math.Pow(roots[i].Abs, n));
                Assert.True(val.Abs < tol * scale,
                    $"Root {i}: |P({roots[i]})| = {val.Abs:G17} exceeds {tol * scale:G17}");
            }
        }

        private static void AssertRootsMatch(Cmplx[] expected, Cmplx* actual, int n, double tol = Tol)
        {
            bool[] matched = new bool[n];
            for (int i = 0; i < n; i++)
            {
                bool found = false;
                for (int j = 0; j < n; j++)
                {
                    if (matched[j]) continue;
                    if ((expected[i] - actual[j]).Abs < tol * Math.Max(1.0, expected[i].Abs))
                    {
                        matched[j] = true;
                        found      = true;
                        break;
                    }
                }

                Assert.True(found,
                    $"Expected root {i} ({expected[i]}) not matched in actual roots");
            }
        }

        // --- degree 5: simplest case that exercises Newton + deflation + quartic ---

        [Fact]
        public void Degree5_RealRoots()
        {
            const int N      = 5;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new Cmplx[] { 1.0, -2.0, 3.0, -4.0, 5.0 };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
            AssertRootsMatch(expected, roots, N);
        }

        [Fact]
        public void Degree5_MixedRealAndComplex()
        {
            const int N      = 5;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new[] { new Cmplx(1.0, 2.0), new Cmplx(1.0, -2.0), new Cmplx(-1.0, 0.5), new Cmplx(-1.0, -0.5), new Cmplx(3.0) };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
            AssertRootsMatch(expected, roots, N);
        }

        // --- degree 6: two Newton + deflate iterations before quartic ---

        [Fact]
        public void Degree6_ConjugatePairs()
        {
            const int N      = 6;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new[] { new Cmplx(2.0, 1.0), new Cmplx(2.0, -1.0), new Cmplx(-1.0, 3.0), new Cmplx(-1.0, -3.0), new Cmplx(0.5, 0.5), new Cmplx(0.5, -0.5) };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
            AssertRootsMatch(expected, roots, N);
        }

        // --- degree 8: unit circle roots, trigonometric mode ---

        [Fact]
        public void Degree8_UnitCircle_Trigonometric()
        {
            const int N      = 8;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new Cmplx[N];
            for (int i = 0; i < N; i++)
                expected[i] = Cmplx.Polar(1.0, 2.0 * Math.PI * i / N + 0.1 * i);
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, true, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
        }

        // --- degree 16: MOID-like trigonometric polynomial ---

        [Fact]
        public void Degree16_TrigonometricSymmetry()
        {
            const int N      = 16;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new Cmplx[N];
            for (int i = 0; i < 8; i++)
            {
                double angle = Math.PI * i / 8.0 + 0.15;
                double r     = 1.0 + 0.02 * Math.Sin(i);
                expected[2 * i]     = Cmplx.Polar(r, angle);
                expected[2 * i + 1] = Cmplx.Polar(1.0 / r, -angle);
            }

            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, true, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots, 1e-8);
        }

        [Fact]
        public void Degree16_AllOnUnitCircle()
        {
            const int N      = 16;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new Cmplx[N];
            for (int i = 0; i < N; i++)
                expected[i] = Cmplx.Polar(1.0, 2.0 * Math.PI * i / N + 0.03 * i);
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, true, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots, 1e-8);
        }

        // --- pre-seeded initial guesses (as MOID_fast does) ---

        [Fact]
        public void Degree16_PreSeededGuesses()
        {
            const int N      = 16;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            var expected = new Cmplx[N];
            for (int i = 0; i < 8; i++)
            {
                double angle = Math.PI * i / 8.0 + 0.15;
                double r     = 1.0 + 0.02 * Math.Sin(i);
                expected[2 * i]     = Cmplx.Polar(r, angle);
                expected[2 * i + 1] = Cmplx.Polar(1.0 / r, -angle);
            }

            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            for (int i = 0; i < 8; i++)
                roots[i] = expected[i] * Cmplx.Polar(1.0, 0.01);
            for (int i = 8; i < N; i++)
                roots[i] = 0.0;

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, true, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots, 1e-8);
        }

        // --- trigonometric guess seeding: verify 1/z* prediction works ---

        [Fact]
        public void TrigonometricGuessSeeding_FindsConjugateInversePair()
        {
            const int N      = 6;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            var z0    = Cmplx.Polar(1.1, 0.5);
            var z0inv = Cmplx.Polar(1.0 / 1.1, -0.5);

            var expected = new[] { z0, z0inv, Cmplx.Polar(1.0, 1.5), Cmplx.Polar(1.0, -1.5), Cmplx.Polar(1.0, 2.5), Cmplx.Polar(1.0, -2.5) };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            roots[0] = z0 * Cmplx.Polar(1.0, 0.01);
            for (int i = 1; i < N; i++) roots[i] = 0.0;

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, true, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
        }

        // --- real-coefficient (non-trigonometric) path ---

        [Fact]
        public void NonTrigonometric_RealCoefficients_ConjugateGuessing()
        {
            const int N      = 6;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new[] { new Cmplx(2.0, 3.0), new Cmplx(2.0, -3.0), new Cmplx(-1.0, 0.5), new Cmplx(-1.0, -0.5), new Cmplx(4.0), new Cmplx(-2.0) };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
            AssertRootsMatch(expected, roots, N);
        }

        [Fact]
        public void NonTrigonometric_NearlyRealRoot_SwapsComponents()
        {
            const int N      = 6;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new[] { new Cmplx(5.0, 1e-5), new Cmplx(5.0, -1e-5), new Cmplx(1.0, 2.0), new Cmplx(1.0, -2.0), new Cmplx(-3.0), new Cmplx(2.0) };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots, 1e-6);
        }

        // --- non-monic polynomial ---

        [Fact]
        public void NonMonicPolynomial()
        {
            const int N      = 6;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new[] { new Cmplx(1.0, 1.0), new Cmplx(1.0, -1.0), new Cmplx(-2.0, 0.5), new Cmplx(-2.0, -0.5), new Cmplx(3.0), new Cmplx(-1.0) };
            BuildPolynomial(expected, c);
            var scale                         = new Cmplx(5.0, -2.0);
            for (int i = 0; i <= N; i++) c[i] = c[i] * scale;
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
            AssertRootsMatch(expected, roots, N);
        }

        // --- exactly degree 5: one Newton step + quartic ---

        [Fact]
        public void Degree5_ExactlyOneDeflation()
        {
            const int N      = 5;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new[] { Cmplx.Polar(1.0, 0.3), Cmplx.Polar(1.0, 1.2), Cmplx.Polar(1.0, 2.5), Cmplx.Polar(1.0, 3.8), Cmplx.Polar(1.0, 5.1) };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, true, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
        }

        // --- iteration count is reasonable ---

        [Fact]
        public void IterationCount_Degree16()
        {
            const int N      = 16;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new Cmplx[N];
            for (int i = 0; i < N; i++)
                expected[i] = Cmplx.Polar(1.0 + 0.01 * Math.Sin(i),
                    2.0 * Math.PI * i / N + 0.03 * i);
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            ulong iters = PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, true, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots, 1e-8);
            Assert.True(iters < 1000, $"Total iterations {iters} seems excessive for degree 16");
        }

        // --- near-double roots (ill-conditioned but should still work) ---

        [Fact]
        public void NearDoubleRoots()
        {
            const int N      = 6;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new[] { new Cmplx(1.0), new Cmplx(1.0 + 1e-4), new Cmplx(-2.0), new Cmplx(-2.0 + 1e-4), new Cmplx(0.5, 1.0), new Cmplx(0.5, -1.0) };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots, 1e-4);
        }

        // --- verify residuals against ORIGINAL (undeflated) coefficients ---

        [Fact]
        public void Degree16_ResidualsAgainstOriginalCoefficients()
        {
            const int N      = 16;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new Cmplx[N];
            for (int i = 0; i < 8; i++)
            {
                double angle = 0.4 + Math.PI * i / 4.0;
                double r     = 1.0 + 0.015 * (i - 4);
                expected[2 * i]     = Cmplx.Polar(r, angle);
                expected[2 * i + 1] = Cmplx.Polar(1.0 / r, -angle);
            }

            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, true, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots, 1e-6);
        }

        // --- mixed magnitude roots ---

        [Fact]
        public void Degree8_MixedMagnitudes()
        {
            const int N      = 8;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            var expected = new[] { new Cmplx(0.1, 0.2), new Cmplx(0.1, -0.2), new Cmplx(0.5), new Cmplx(-0.8), new Cmplx(2.0, 1.0), new Cmplx(2.0, -1.0), new Cmplx(5.0), new Cmplx(-3.0) };
            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots);
            AssertRootsMatch(expected, roots, N);
        }

        // --- stress deflation with many roots ---

        [Fact]
        public void Degree12_StressDeflation()
        {
            const int N      = 12;
            Cmplx*    c      = stackalloc Cmplx[N + 1];
            Cmplx*    cSaved = stackalloc Cmplx[N + 1];
            Cmplx*    roots  = stackalloc Cmplx[N];

            for (int i = 0; i < N; i++) roots[i] = 0.0;

            // 6 conjugate pairs — gives real coefficients, matching trigonometric=false
            var expected = new Cmplx[N];
            for (int i = 0; i < 6; i++)
            {
                double r     = 0.8 + 0.4 * (i % 3);
                double angle = 2.0 * Math.PI * i / 6 + 0.05;
                expected[2 * i]     = Cmplx.Polar(r, angle);
                expected[2 * i + 1] = Cmplx.Polar(r, -angle);
            }

            BuildPolynomial(expected, c);
            CopyCoeffs(c, cSaved, N);

            PolySolver.PolynomialRoots(N, c, roots, 1e-14, 1e-8, false, PolySolver.MakeRng());

            AssertAllResiduals(cSaved, N, roots, 1e-6);
        }
    }
}
