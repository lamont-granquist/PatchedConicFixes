using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    namespace PatchedConicFixes.Tests
    {
        public unsafe class SolveQuarticTests
        {
            // Build monic quartic from four roots: (z - z0)(z - z1)(z - z2)(z - z3)
            private static void BuildQuartic(Cmplx z0, Cmplx z1, Cmplx z2, Cmplx z3, Cmplx* c)
            {
                // expand via symmetric polynomials
                Cmplx e1 = z0 + z1 + z2 + z3;
                Cmplx e2 = z0 * z1 + z0 * z2 + z0 * z3 + z1 * z2 + z1 * z3 + z2 * z3;
                Cmplx e3 = z0 * z1 * z2 + z0 * z1 * z3 + z0 * z2 * z3 + z1 * z2 * z3;
                Cmplx e4 = z0 * z1 * z2 * z3;

                c[4] = 1.0;
                c[3] = -e1;
                c[2] = e2;
                c[1] = -e3;
                c[0] = e4;
            }

            private static Cmplx Eval(Cmplx* c, Cmplx z)
                => (((c[4] * z + c[3]) * z + c[2]) * z + c[1]) * z + c[0];

            private static void AssertResidualSmall(Cmplx* c, Cmplx root, double tol = 1e-10)
            {
                Cmplx  val   = Eval(c, root);
                double scale = Math.Max(1.0, c[4].Abs * root.Norm * root.Norm);
                Assert.True(val.Abs < tol * scale,
                    $"Residual |P({root.Re:G17} + {root.Im:G17}i)| = {val.Abs:G17} exceeds {tol * scale:G17}");
            }

            private static void AssertAllResiduals(Cmplx* c, Cmplx* roots, double tol = 1e-10)
            {
                for (int i = 0; i < 4; i++)
                    AssertResidualSmall(c, roots[i], tol);
            }

            // Check that two sets of 4 roots match in some order
            private static void AssertRootsMatch(Cmplx[] expected, Cmplx* actual, double tol = 1e-10)
            {
                bool[] matched = new bool[4];
                for (int i = 0; i < 4; i++)
                {
                    bool found = false;
                    for (int j = 0; j < 4; j++)
                    {
                        if (matched[j]) continue;
                        if (Close(expected[i], actual[j], tol))
                        {
                            matched[j] = true;
                            found      = true;
                            break;
                        }
                    }

                    Assert.True(found,
                        $"Expected root [{expected[i].Re:G17}, {expected[i].Im:G17}] not found in actual roots");
                }
            }

            private static bool Close(Cmplx a, Cmplx b, double tol)
                => Math.Abs(a.Re - b.Re) < tol && Math.Abs(a.Im - b.Im) < tol;

            [Fact]
            public void FourDistinctRealRoots()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                var z = new Cmplx[] { 1.0, -2.0, 3.0, -4.0 };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots);
                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void TwoConjugatePairs()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z = { new Cmplx(1.0, 2.0), new Cmplx(1.0, -2.0), new Cmplx(-3.0, 1.0), new Cmplx(-3.0, -1.0) };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots);
                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void OneConjugatePairTwoReal()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z = { new Cmplx(0.5, 1.5), new Cmplx(0.5, -1.5), new Cmplx(2.0), new Cmplx(-1.0) };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots);
                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void RootsOnUnitCircle()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                // Four evenly spaced on |z|=1: e^(ikπ/4) for k = 1,3,5,7
                Cmplx[] z = { Cmplx.Polar(1.0, Math.PI / 4.0), Cmplx.Polar(1.0, 3.0 * Math.PI / 4.0), Cmplx.Polar(1.0, 5.0 * Math.PI / 4.0), Cmplx.Polar(1.0, 7.0 * Math.PI / 4.0) };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots);
                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void RootsNearUnitCircle()
            {
                // Simulates the MOID case: roots slightly off |z|=1
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z =
                {
                    Cmplx.Polar(1.001, 0.3), Cmplx.Polar(0.999, 0.3), // near-inverse of above
                    Cmplx.Polar(1.002, 2.1), Cmplx.Polar(0.998, 2.1)
                };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots);
                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void DoubleRealRoot()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                // (z-2)²(z-5)(z+1)
                BuildQuartic(2.0, 2.0, 5.0, -1.0, c);

                PolySolver.SolveQuartic(c, roots);

                // Relaxed tolerance — double roots are ill-conditioned
                AssertAllResiduals(c, roots, 1e-6);
            }

            [Fact]
            public void TwoDoubleRoots()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                // (z-1)²(z+3)²
                BuildQuartic(1.0, 1.0, -3.0, -3.0, c);

                PolySolver.SolveQuartic(c, roots);

                AssertAllResiduals(c, roots, 1e-6);
            }

            [Fact]
            public void TripleRoot()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                // (z-4)³(z+2)
                BuildQuartic(4.0, 4.0, 4.0, -2.0, c);

                PolySolver.SolveQuartic(c, roots);

                AssertAllResiduals(c, roots, 1e-4);
            }

            [Fact]
            public void QuadrupleRoot()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                // (z-3)⁴
                BuildQuartic(3.0, 3.0, 3.0, 3.0, c);

                PolySolver.SolveQuartic(c, roots);

                AssertAllResiduals(c, roots, 1e-3);
            }

            [Fact]
            public void PureImaginaryRoots()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z = { new Cmplx(0.0, 1.0), new Cmplx(0.0, -1.0), new Cmplx(0.0, 3.0), new Cmplx(0.0, -3.0) };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots);
                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void MixedInsideOutsideUnitCircle()
            {
                // Two roots inside, two outside — exercises both polish branches
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z =
                {
                    new Cmplx(0.1, 0.2),  // |z| ≈ 0.22
                    new Cmplx(-0.3, 0.4), // |z| = 0.5
                    new Cmplx(3.0, -4.0), // |z| = 5
                    new Cmplx(-2.0, 1.0)  // |z| ≈ 2.24
                };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots);
                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void VerySmallRoots()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z = { new Cmplx(1e-4, 2e-4), new Cmplx(1e-4, -2e-4), new Cmplx(-3e-4, 1e-4), new Cmplx(-3e-4, -1e-4) };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots, 1e-8);
                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void VeryLargeRoots()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z = { new Cmplx(1e4, 2e4), new Cmplx(1e4, -2e4), new Cmplx(-3e4), new Cmplx(5e4) };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertAllResiduals(c, roots, 1e-6);
            }

            [Fact]
            public void NonMonicPolynomial()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z     = { new Cmplx(1.0, 1.0), new Cmplx(1.0, -1.0), new Cmplx(-2.0), new Cmplx(0.5) };
                var     scale = new Cmplx(3.0, -7.0);

                // Build monic first, then multiply all coefficients by scale
                BuildQuartic(z[0], z[1], z[2], z[3], c);
                for (int i = 0; i < 5; i++) c[i] = c[i] * scale;

                PolySolver.SolveQuartic(c, roots);

                AssertRootsMatch(z, roots);
                AssertAllResiduals(c, roots);
            }

            // NOTE: Ferrari's method is poor for widely separated magnitudes.  For a general purpose
            // polynomial root finding you would either want to reduce the polynomial with 2 additional
            // Newton solves and use the quadratic solver -- or at least check the result from
            // Ferrari's method and if the residuals were poor then to fall back to quadratic solving.  For
            // a MOID computation the roots should be near the unit circle and therefore it shouldn't matter.
            [Fact]
            public void WidelySeparatedMagnitudes()
            {
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                var z = new Cmplx[]
                {
                    new Cmplx(1e-3),
                    new Cmplx(-1e-3),
                    new Cmplx(1e3),
                    new Cmplx(-1e3),
                };
                BuildQuartic(z[0], z[1], z[2], z[3], c);

                PolySolver.SolveQuartic(c, roots);

                AssertAllResiduals(c, roots, 1e-6);
            }

            [Fact]
            public void ComplexCoefficients()
            {
                // Polynomial with fully complex coefficients (no conjugation symmetry)
                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                Cmplx[] z = { new Cmplx(1.0, 2.0), new Cmplx(-1.0, 3.0), new Cmplx(2.0, -1.0), new Cmplx(-3.0, -2.0) };

                // Scale by a complex leading coefficient
                var lead = new Cmplx(2.0, 5.0);
                BuildQuartic(z[0], z[1], z[2], z[3], c);
                for (int i = 0; i < 5; i++) c[i] = c[i] * lead;

                PolySolver.SolveQuartic(c, roots);

                AssertAllResiduals(c, roots);
            }

            [Fact]
            public void ResidualSweep()
            {
                // Several cases, just verify P(root) ≈ 0 for all
                Cmplx[][] cases = { new Cmplx[] { 1.0, 2.0, 3.0, 4.0 }, new Cmplx[] { -1.0, -2.0, -3.0, -4.0 }, new[] { new Cmplx(0.5, 0.5), new Cmplx(0.5, -0.5), new Cmplx(-0.5, 0.5), new Cmplx(-0.5, -0.5) }, new[] { Cmplx.Polar(1.0, 0.1), Cmplx.Polar(1.0, 1.2), Cmplx.Polar(1.0, 2.5), Cmplx.Polar(1.0, 4.0) }, new Cmplx[] { 0.0, 0.0, 0.0, 5.0 } };

                Cmplx* c     = stackalloc Cmplx[5];
                Cmplx* roots = stackalloc Cmplx[4];

                foreach (Cmplx[] z in cases)
                {
                    BuildQuartic(z[0], z[1], z[2], z[3], c);
                    PolySolver.SolveQuartic(c, roots);
                    AssertAllResiduals(c, roots, 1e-6);
                }
            }
        }
    }
}
