using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public unsafe class SolveQuadraticTests
    {
        // Build c[0] + c[1]z + c[2]z² from roots z0, z1 with c[2] = 1.
        private static void BuildQuadratic(Cmplx z0, Cmplx z1, Cmplx* c)
        {
            c[2] = new Cmplx(1.0);
            c[1] = -(z0 + z1);
            c[0] = z0 * z1;
        }

        // Evaluate c[0] + c[1]z + c[2]z²
        private static Cmplx Eval(Cmplx* c, Cmplx z) => (c[2] * z + c[1]) * z + c[0];

        private static void AssertRootsMatch(Cmplx expected0, Cmplx expected1,
            Cmplx actual0, Cmplx actual1,
            double tol = 1e-12)
        {
            bool order1 = Close(expected0, actual0, tol) && Close(expected1, actual1, tol);
            bool order2 = Close(expected0, actual1, tol) && Close(expected1, actual0, tol);
            Assert.True(order1 || order2,
                $"Roots [{actual0.Re:G17}, {actual0.Im:G17}] and [{actual1.Re:G17}, {actual1.Im:G17}] " +
                $"don't match expected [{expected0.Re:G17}, {expected0.Im:G17}] and [{expected1.Re:G17}, {expected1.Im:G17}]");
        }

        private static bool Close(Cmplx a, Cmplx b, double tol)
            => Math.Abs(a.Re - b.Re) < tol && Math.Abs(a.Im - b.Im) < tol;

        private static void AssertResidualSmall(Cmplx* c, Cmplx root, double tol = 1e-12)
        {
            Cmplx  val   = Eval(c, root);
            double scale = Math.Max(1.0, c[2].Abs * root.Norm);
            Assert.True(val.Abs < tol * scale,
                $"Residual |P({root.Re:G17} + {root.Im:G17}i)| = {val.Abs:G17} exceeds {tol * scale:G17}");
        }

        [Fact]
        public void TwoDistinctRealRoots()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0 = new Cmplx(2.0);
            var z1 = new Cmplx(-3.0);
            BuildQuadratic(z0, z1, c);

            PolySolver.SolveQuadratic(c, roots);

            AssertRootsMatch(z0, z1, roots[0], roots[1]);
        }

        [Fact]
        public void ComplexConjugateRoots()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0 = new Cmplx(1.0, 3.0);
            var z1 = new Cmplx(1.0, -3.0);
            BuildQuadratic(z0, z1, c);

            PolySolver.SolveQuadratic(c, roots);

            AssertRootsMatch(z0, z1, roots[0], roots[1]);
        }

        [Fact]
        public void DoubleRoot()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0 = new Cmplx(5.0);
            BuildQuadratic(z0, z0, c);

            PolySolver.SolveQuadratic(c, roots);

            // Both roots should be near the double root
            AssertRootsMatch(z0, z0, roots[0], roots[1], 1e-8);
        }

        [Fact]
        public void RootsOnUnitCircle()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            // e^(i·π/3) and e^(-i·π/3) — conjugate pair on |z|=1
            var z0 = Cmplx.Polar(1.0, Math.PI / 3.0);
            var z1 = Cmplx.Polar(1.0, -Math.PI / 3.0);
            BuildQuadratic(z0, z1, c);

            PolySolver.SolveQuadratic(c, roots);

            AssertRootsMatch(z0, z1, roots[0], roots[1]);
        }

        [Fact]
        public void OneRootInsideOneOutsideUnitCircle()
        {
            // Exercises both branches of the polish step
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0 = new Cmplx(0.1, 0.2);  // |z| ≈ 0.22
            var z1 = new Cmplx(3.0, -4.0); // |z| = 5
            BuildQuadratic(z0, z1, c);

            PolySolver.SolveQuadratic(c, roots);

            AssertRootsMatch(z0, z1, roots[0], roots[1]);
        }

        [Fact]
        public void VerySmallRoots()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0 = new Cmplx(1e-8, 1e-8);
            var z1 = new Cmplx(-1e-8, 2e-8);
            BuildQuadratic(z0, z1, c);

            PolySolver.SolveQuadratic(c, roots);

            AssertRootsMatch(z0, z1, roots[0], roots[1], 1e-10);
        }

        [Fact]
        public void VeryLargeRoots()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0 = new Cmplx(1e8, -1e8);
            var z1 = new Cmplx(-1e8, 1e7);
            BuildQuadratic(z0, z1, c);

            PolySolver.SolveQuadratic(c, roots);

            AssertRootsMatch(z0, z1, roots[0], roots[1], 1e-4);
        }

        [Fact]
        public void NonMonicPolynomial()
        {
            // Leading coefficient ≠ 1
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0    = new Cmplx(1.0, 2.0);
            var z1    = new Cmplx(-1.0, 3.0);
            var scale = new Cmplx(7.0, -3.0);

            // scale·(z - z0)(z - z1)
            c[2] = scale;
            c[1] = scale * -(z0 + z1);
            c[0] = scale * (z0 * z1);

            PolySolver.SolveQuadratic(c, roots);

            AssertRootsMatch(z0, z1, roots[0], roots[1]);
        }

        [Fact]
        public void PureImaginaryRoots()
        {
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0 = new Cmplx(0.0, 4.0);
            var z1 = new Cmplx(0.0, -4.0);
            BuildQuadratic(z0, z1, c);

            PolySolver.SolveQuadratic(c, roots);

            AssertRootsMatch(z0, z1, roots[0], roots[1]);
        }

        [Fact]
        public void WidelySeperatedMagnitudes()
        {
            // One root near zero, one far away — tests Vieta's formula path
            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            var z0 = new Cmplx(1e-15);
            var z1 = new Cmplx(1e15);
            BuildQuadratic(z0, z1, c);

            PolySolver.SolveQuadratic(c, roots);

            // Residuals should be small even if root matching is hard
            AssertResidualSmall(c, roots[0]);
            AssertResidualSmall(c, roots[1]);
        }

        [Fact]
        public void ResidualCheckOnAllRoots()
        {
            // Fuzz-like: several random-ish polynomials, verify P(root) ≈ 0
            var cases = new (double r0, double i0, double r1, double i1)[] { (1.0, 0.0, -1.0, 0.0), (0.5, 0.5, 0.5, -0.5), (2.3, -1.7, -0.8, 4.1), (0.001, 0.0, 1000.0, 0.0), (0.0, 1.0, 0.0, -1.0) };

            Cmplx* c     = stackalloc Cmplx[3];
            Cmplx* roots = stackalloc Cmplx[2];

            foreach ((double r0, double i0, double r1, double i1) in cases)
            {
                BuildQuadratic(new Cmplx(r0, i0), new Cmplx(r1, i1), c);
                PolySolver.SolveQuadratic(c, roots);
                AssertResidualSmall(c, roots[0]);
                AssertResidualSmall(c, roots[1]);
            }
        }
    }
}
