using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    namespace PatchedConicFixes.Tests
    {
        public unsafe class ExtractLinearFactorTests
        {
            // Build monic polynomial from roots: (z - z0)(z - z1)...(z - zn)
            // Coefficients in ascending order: c[0] + c[1]z + ... + c[n]z^n
            private static void BuildPolynomial(Cmplx[] roots, Cmplx* c)
            {
                int n = roots.Length;
                // start with c[0] = 1 (constant polynomial "1")
                for (int i = 0; i <= n; i++) c[i] = 0.0;
                c[0] = 1.0;

                // multiply by (z - root) one at a time
                for (int k = 0; k < n; k++)
                {
                    // expand c * (z - roots[k]) in place, top-down
                    c[k + 1] = c[k];
                    for (int i = k; i > 0; i--)
                        c[i] = c[i - 1] - roots[k] * c[i];
                    c[0] = -roots[k] * c[0];
                }
            }

            private static Cmplx EvalPoly(Cmplx* c, int degree, Cmplx z)
            {
                Cmplx val = c[degree];
                for (int i = degree - 1; i >= 0; i--)
                    val = val * z + c[i];
                return val;
            }

            private static void AssertResidualSmall(Cmplx* c, int degree, Cmplx root, double tol = 1e-10)
            {
                Cmplx  val   = EvalPoly(c, degree, root);
                double scale = Math.Max(1.0, root.Abs);
                Assert.True(val.Abs < tol * Math.Pow(scale, degree),
                    $"|P({root.Re:G17} + {root.Im:G17}i)| = {val.Abs:G17}");
            }

            [Fact]
            public void ExtractRootInsideUnitCircle_PointerAdvances()
            {
                // (z - 0.5)(z - 3) = z² - 3.5z + 1.5
                Cmplx*  c     = stackalloc Cmplx[3];
                Cmplx[] roots = { new Cmplx(0.5), new Cmplx(3.0) };
                BuildPolynomial(roots, c);

                Cmplx* c_ = PolySolver.ExtractLinearFactor(2, c, new Cmplx(0.5));

                // |root| < 1 so pointer should advance
                Assert.True(c_ == c + 1);

                // Deflated polynomial (degree 1) should vanish at the remaining root
                AssertResidualSmall(c_, 1, new Cmplx(3.0));
            }

            [Fact]
            public void ExtractRootOutsideUnitCircle_PointerStays()
            {
                // (z - 0.5)(z - 3) = z² - 3.5z + 1.5
                Cmplx*  c     = stackalloc Cmplx[3];
                Cmplx[] roots = { new Cmplx(0.5), new Cmplx(3.0) };
                BuildPolynomial(roots, c);

                Cmplx* c_ = PolySolver.ExtractLinearFactor(2, c, new Cmplx(3.0));

                // |root| > 1 so pointer should stay
                Assert.True(c_ == c);

                // Deflated polynomial (degree 1) should vanish at the remaining root
                AssertResidualSmall(c_, 1, new Cmplx(0.5));
            }

            [Fact]
            public void ExtractRootOnUnitCircle()
            {
                Cmplx*  c     = stackalloc Cmplx[3];
                var     z0    = Cmplx.Polar(1.0, Math.PI / 3.0);
                var     z1    = new Cmplx(2.0, -1.0);
                Cmplx[] roots = { z0, z1 };
                BuildPolynomial(roots, c);

                // |z0| = 1.0 exactly, so rootNorm <= 1.0 branch, pointer advances
                Cmplx* c_ = PolySolver.ExtractLinearFactor(2, c, z0);

                Assert.True(c_ == c + 1);
                AssertResidualSmall(c_, 1, z1);
            }

            [Fact]
            public void ExtractRootAtZero()
            {
                // z(z - 2)(z - 3) = z³ - 5z² + 6z
                Cmplx* c     = stackalloc Cmplx[4];
                var    roots = new Cmplx[] { 0.0, 2.0, 3.0 };
                BuildPolynomial(roots, c);

                Cmplx* c_ = PolySolver.ExtractLinearFactor(3, c, 0.0);

                Assert.True(c_ == c + 1);
                AssertResidualSmall(c_, 2, new Cmplx(2.0));
                AssertResidualSmall(c_, 2, new Cmplx(3.0));
            }

            [Fact]
            public void ExtractComplexRootInsideUnitCircle()
            {
                Cmplx*  c     = stackalloc Cmplx[3];
                var     z0    = new Cmplx(0.3, 0.4);  // |z| = 0.5
                var     z1    = new Cmplx(-2.0, 1.0); // |z| ≈ 2.24
                Cmplx[] roots = { z0, z1 };
                BuildPolynomial(roots, c);

                Cmplx* c_ = PolySolver.ExtractLinearFactor(2, c, z0);

                Assert.True(c_ == c + 1);
                AssertResidualSmall(c_, 1, z1);
            }

            [Fact]
            public void ExtractComplexRootOutsideUnitCircle()
            {
                Cmplx*  c     = stackalloc Cmplx[3];
                var     z0    = new Cmplx(0.3, 0.4);  // |z| = 0.5
                var     z1    = new Cmplx(-2.0, 1.0); // |z| ≈ 2.24
                Cmplx[] roots = { z0, z1 };
                BuildPolynomial(roots, c);

                Cmplx* c_ = PolySolver.ExtractLinearFactor(2, c, z1);

                Assert.True(c_ == c);
                AssertResidualSmall(c_, 1, z0);
            }

            [Fact]
            public void SequentialExtraction_Degree4()
            {
                Cmplx* c = stackalloc Cmplx[5];
                Cmplx[] roots =
                {
                    new Cmplx(0.2, 0.3),  // inside
                    new Cmplx(0.2, -0.3), // inside
                    new Cmplx(3.0, 1.0),  // outside
                    new Cmplx(-2.0)       // outside
                };
                BuildPolynomial(roots, c);

                // Extract one at a time, verify remaining roots after each
                Cmplx* c_ = PolySolver.ExtractLinearFactor(4, c, roots[0]);
                AssertResidualSmall(c_, 3, roots[1]);
                AssertResidualSmall(c_, 3, roots[2]);
                AssertResidualSmall(c_, 3, roots[3]);

                c_ = PolySolver.ExtractLinearFactor(3, c_, roots[2]);
                AssertResidualSmall(c_, 2, roots[1]);
                AssertResidualSmall(c_, 2, roots[3]);

                c_ = PolySolver.ExtractLinearFactor(2, c_, roots[1]);
                AssertResidualSmall(c_, 1, roots[3]);
            }

            [Fact]
            public void SequentialExtraction_AlternatingInsideOutside()
            {
                // Alternate extracting inside/outside roots to stress pointer tracking
                Cmplx* c = stackalloc Cmplx[5];
                Cmplx[] roots =
                {
                    new Cmplx(0.1), // inside, ptr advances
                    new Cmplx(5.0), // outside, ptr stays
                    new Cmplx(0.5), // inside, ptr advances
                    new Cmplx(2.0)  // outside, ptr stays
                };
                BuildPolynomial(roots, c);

                Cmplx* c_ = c;
                for (int i = 0; i < 4; i++)
                {
                    c_ = PolySolver.ExtractLinearFactor(4 - i, c_, roots[i]);
                    for (int j = i + 1; j < 4; j++)
                        AssertResidualSmall(c_, 4 - i - 1, roots[j]);
                }
            }

            [Fact]
            public void SequentialExtraction_Degree16_UnitCircleRoots()
            {
                // Simulates the MOID case: 16 roots near |z| = 1
                const int N     = 16;
                Cmplx*    c     = stackalloc Cmplx[N + 1];
                var       roots = new Cmplx[N];
                for (int i = 0; i < N; i++)
                    roots[i] = Cmplx.Polar(1.0, 2.0 * Math.PI * i / N);
                BuildPolynomial(roots, c);

                Cmplx* c_ = c;
                for (int i = 0; i < N - 2; i++)
                {
                    c_ = PolySolver.ExtractLinearFactor(N - i, c_, roots[i]);

                    // Spot-check a few remaining roots
                    AssertResidualSmall(c_, N - i - 1, roots[N - 1], 1e-6);
                }

                // Final quadratic should have the last two roots
                AssertResidualSmall(c_, 2, roots[N - 2], 1e-4);
                AssertResidualSmall(c_, 2, roots[N - 1], 1e-4);
            }

            [Fact]
            public void DoubleRoot_ExtractOneInstance()
            {
                // (z - 2)²(z - 5)
                Cmplx* c     = stackalloc Cmplx[4];
                var    roots = new Cmplx[] { 2.0, 2.0, 5.0 };
                BuildPolynomial(roots, c);

                Cmplx* c_ = PolySolver.ExtractLinearFactor(3, c, new Cmplx(2.0));

                // Remaining quadratic should still have roots at 2 and 5
                AssertResidualSmall(c_, 2, new Cmplx(2.0));
                AssertResidualSmall(c_, 2, new Cmplx(5.0));
            }

            [Fact]
            public void NonMonicPolynomial()
            {
                // 7·(z - 1)(z + 2)
                Cmplx* c     = stackalloc Cmplx[3];
                var    roots = new Cmplx[] { 1.0, -2.0 };
                BuildPolynomial(roots, c);
                var scale = new Cmplx(7.0, -3.0);
                for (int i = 0; i < 3; i++)
                    c[i] = c[i] * scale;

                Cmplx* c_ = PolySolver.ExtractLinearFactor(2, c, new Cmplx(1.0));

                // Deflated linear poly should vanish at -2
                AssertResidualSmall(c_, 1, new Cmplx(-2.0));
            }
        }
    }
}
