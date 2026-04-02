using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public unsafe class PolynomialValDerTests
    {
        private static void AssertClose(Cmplx expected, Cmplx actual, double tol = 1e-12, string label = "")
        {
            Assert.True(Math.Abs(expected.Re - actual.Re) < tol &&
                Math.Abs(expected.Im - actual.Im) < tol,
                $"{label} Expected ({expected.Re:G17}, {expected.Im:G17}) " +
                $"got ({actual.Re:G17}, {actual.Im:G17})");
        }

        // Evaluate derivative by finite difference for cross-checking
        private static Cmplx FiniteDiffDer(int n, Cmplx* c, Cmplx z, bool forward)
        {
            double h = 1e-7;
            PolySolver.PolynomialValDer(n, c, z + h, out Cmplx fp, out _, forward);
            PolySolver.PolynomialValDer(n, c, z - h, out Cmplx fm, out _, forward);
            return (fp - fm) / (2.0 * h);
        }

        // --- constant polynomial: c[0] = k ---

        [Fact]
        public void Constant_Forward()
        {
            Cmplx* c = stackalloc Cmplx[1];
            c[0] = new Cmplx(7.0, -3.0);

            PolySolver.PolynomialValDer(0, c, new Cmplx(99.0), out Cmplx val, out Cmplx der, true);

            AssertClose(c[0], val, label: "val");
            AssertClose(0.0, der, label: "der");
        }

        [Fact]
        public void Constant_Backward()
        {
            Cmplx* c = stackalloc Cmplx[1];
            c[0] = new Cmplx(7.0, -3.0);

            PolySolver.PolynomialValDer(0, c, new Cmplx(99.0), out Cmplx val, out Cmplx der, false);

            AssertClose(c[0], val, label: "val");
            AssertClose(0.0, der, label: "der");
        }

        // --- linear polynomial: c[0] + c[1]z ---

        [Fact]
        public void Linear_Forward()
        {
            // P(z) = 3 + 5z, P'(z) = 5
            Cmplx* c = stackalloc Cmplx[2];
            c[0] = 3.0;
            c[1] = 5.0;
            Cmplx z = 2.0;

            PolySolver.PolynomialValDer(1, c, z, out Cmplx val, out Cmplx der, true);

            AssertClose(13.0, val, label: "val");
            AssertClose(5.0, der, label: "der");
        }

        [Fact]
        public void Linear_Backward()
        {
            // Q(w) = c[1] + c[0]w = 5 + 3w, Q'(w) = 3
            Cmplx* c = stackalloc Cmplx[2];
            c[0] = 3.0;
            c[1] = 5.0;
            Cmplx w = 2.0;

            PolySolver.PolynomialValDer(1, c, w, out Cmplx val, out Cmplx der, false);

            AssertClose(11.0, val, label: "val");
            AssertClose(3.0, der, label: "der");
        }

        // --- quadratic: c[0] + c[1]z + c[2]z² ---

        [Fact]
        public void Quadratic_Forward()
        {
            // P(z) = 1 + 2z + 3z², P'(z) = 2 + 6z
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = 1.0;
            c[1] = 2.0;
            c[2] = 3.0;
            Cmplx z = 4.0;

            PolySolver.PolynomialValDer(2, c, z, out Cmplx val, out Cmplx der, true);

            // P(4) = 1 + 8 + 48 = 57
            AssertClose(57.0, val, label: "val");
            // P'(4) = 2 + 24 = 26
            AssertClose(26.0, der, label: "der");
        }

        [Fact]
        public void Quadratic_Backward()
        {
            // Q(w) = c[2] + c[1]w + c[0]w² = 3 + 2w + 1w²
            // Q'(w) = 2 + 2w
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = 1.0;
            c[1] = 2.0;
            c[2] = 3.0;
            Cmplx w = 4.0;

            PolySolver.PolynomialValDer(2, c, w, out Cmplx val, out Cmplx der, false);

            // Q(4) = 3 + 8 + 16 = 27
            AssertClose(27.0, val, label: "val");
            // Q'(4) = 2 + 8 = 10
            AssertClose(10.0, der, label: "der");
        }

        // --- cubic: known values ---

        [Fact]
        public void Cubic_Forward()
        {
            // P(z) = 2 - z + 3z² + z³
            // P'(z) = -1 + 6z + 3z²
            Cmplx* c = stackalloc Cmplx[4];
            c[0] = 2.0;
            c[1] = -1.0;
            c[2] = 3.0;
            c[3] = 1.0;
            Cmplx z = 2.0;

            PolySolver.PolynomialValDer(3, c, z, out Cmplx val, out Cmplx der, true);

            // P(2) = 2 - 2 + 12 + 8 = 20
            AssertClose(20.0, val, label: "val");
            // P'(2) = -1 + 12 + 12 = 23
            AssertClose(23.0, der, label: "der");
        }

        // --- complex coefficients and argument ---

        [Fact]
        public void ComplexCoefficients_Forward()
        {
            // P(z) = (1+i) + (2-i)z + (3+2i)z²
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = new Cmplx(1.0, 1.0);
            c[1] = new Cmplx(2.0, -1.0);
            c[2] = new Cmplx(3.0, 2.0);
            var z = new Cmplx(1.0, 1.0);

            PolySolver.PolynomialValDer(2, c, z, out Cmplx val, out Cmplx der, true);

            // z² = (1+i)² = 2i
            // P(z) = (1+i) + (2-i)(1+i) + (3+2i)(2i)
            //       = (1+i) + (3+i) + (-4+6i)
            //       = (0 + 8i)
            AssertClose(new Cmplx(0.0, 8.0), val, label: "val");

            // P'(z) = (2-i) + 2(3+2i)(1+i)
            //        = (2-i) + 2(1+8i) ... wait let me recompute
            // 2*(3+2i)*(1+i) = 2*(3+3i+2i-2) = 2*(1+5i) = (2+10i)
            // P'(z) = (2-i) + (2+10i) = (4+9i)
            AssertClose(new Cmplx(4.0, 9.0), der, label: "der");
        }

        [Fact]
        public void ComplexCoefficients_Backward()
        {
            // Same polynomial, backward evaluation at w
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = new Cmplx(1.0, 1.0);
            c[1] = new Cmplx(2.0, -1.0);
            c[2] = new Cmplx(3.0, 2.0);
            var w = new Cmplx(0.5, -0.5);

            PolySolver.PolynomialValDer(2, c, w, out Cmplx val, out Cmplx der, false);

            // Cross-check with finite differences
            Cmplx derFD = FiniteDiffDer(2, c, w, false);

            // Recompute directly: Q(w) = c[0]w² + c[1]w + c[2]
            Cmplx w2          = w * w;
            Cmplx expectedVal = c[0] * w2 + c[1] * w + c[2];
            // Q'(w) = 2*c[0]*w + c[1]
            Cmplx expectedDer = 2.0 * c[0] * w + c[1];

            AssertClose(expectedVal, val, label: "val");
            AssertClose(expectedDer, der, label: "der");
            AssertClose(expectedDer, derFD, 1e-5, "der vs FD");
        }

        // --- forward and backward agree at z=1 for real coefficients ---

        [Fact]
        public void ForwardBackward_AgreeAtZEquals1()
        {
            // When z = w = 1, forward and backward evaluate the same sums
            Cmplx* c = stackalloc Cmplx[5];
            c[0] = 2.0;
            c[1] = -3.0;
            c[2] = 1.0;
            c[3] = 4.0;
            c[4] = -1.0;

            PolySolver.PolynomialValDer(4, c, 1.0, out Cmplx valF, out Cmplx derF, true);
            PolySolver.PolynomialValDer(4, c, 1.0, out Cmplx valB, out Cmplx derB, false);

            // Values should be identical
            AssertClose(valF, valB, label: "val");
            // Derivatives differ: forward is P'(1), backward is Q'(1)
            // But for real symmetric c: they happen to agree when c is palindromic.
            // In general they differ, so just cross-check each against known values.

            // P(1) = 2 - 3 + 1 + 4 - 1 = 3
            AssertClose(3.0, valF, label: "valF");
            // P'(1) = -3 + 2 + 12 - 4 = 7
            AssertClose(7.0, derF, label: "derF");
            // Q(1) = c[4] + c[3] + c[2] + c[1] + c[0] = 3 (same sum)
            AssertClose(3.0, valB, label: "valB");
            // Q'(1) = 4*c[4] + 3*c[3] + 2*c[2] + c[1] ... wait, backward:
            // Q'(w) = 4*c[0]*w³ + 3*c[1]*w² + 2*c[2]*w + c[3]
            // Q'(1) = 4*2 + 3*(-3) + 2*1 + 4 = 8 - 9 + 2 + 4 = 5
            AssertClose(5.0, derB, label: "derB");
        }

        // --- derivative matches finite difference ---

        [Fact]
        public void DerivativeMatchesFiniteDifference_Forward()
        {
            Cmplx* c = stackalloc Cmplx[6];
            c[0] = new Cmplx(1.0, -2.0);
            c[1] = new Cmplx(3.0, 1.0);
            c[2] = new Cmplx(-1.0, 4.0);
            c[3] = new Cmplx(2.0, -3.0);
            c[4] = new Cmplx(0.5, 0.5);
            c[5] = new Cmplx(-1.0, 1.0);
            var z = new Cmplx(0.7, -0.3);

            PolySolver.PolynomialValDer(5, c, z, out _, out Cmplx der, true);
            Cmplx derFD = FiniteDiffDer(5, c, z, true);

            AssertClose(der, derFD, 1e-5, "forward der vs FD");
        }

        [Fact]
        public void DerivativeMatchesFiniteDifference_Backward()
        {
            Cmplx* c = stackalloc Cmplx[6];
            c[0] = new Cmplx(1.0, -2.0);
            c[1] = new Cmplx(3.0, 1.0);
            c[2] = new Cmplx(-1.0, 4.0);
            c[3] = new Cmplx(2.0, -3.0);
            c[4] = new Cmplx(0.5, 0.5);
            c[5] = new Cmplx(-1.0, 1.0);
            var z = new Cmplx(0.7, -0.3);

            PolySolver.PolynomialValDer(5, c, z, out _, out Cmplx der, false);
            Cmplx derFD = FiniteDiffDer(5, c, z, false);

            AssertClose(der, derFD, 1e-5, "backward der vs FD");
        }

        // --- evaluating at a known root gives val ≈ 0, der ≠ 0 ---

        [Fact]
        public void AtKnownRoot_ValIsZero_DerIsNonzero()
        {
            // (z - 2)(z - 3) = z² - 5z + 6
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = 6.0;
            c[1] = -5.0;
            c[2] = 1.0;

            PolySolver.PolynomialValDer(2, c, 2.0, out Cmplx val, out Cmplx der, true);

            AssertClose(0.0, val, label: "val at root");
            // P'(2) = 2*2 - 5 = -1
            AssertClose(-1.0, der, label: "der at root");
        }

        // --- at a double root, both val and der are zero ---

        [Fact]
        public void AtDoubleRoot_BothZero()
        {
            // (z - 2)² = z² - 4z + 4
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = 4.0;
            c[1] = -4.0;
            c[2] = 1.0;

            PolySolver.PolynomialValDer(2, c, 2.0, out Cmplx val, out Cmplx der, true);

            AssertClose(0.0, val, label: "val at double root");
            AssertClose(0.0, der, label: "der at double root");
        }

        // --- degree 16 with unit circle roots (MOID-like) ---

        [Fact]
        public void Degree16_FiniteDifferenceCheck()
        {
            const int N = 16;
            Cmplx*    c = stackalloc Cmplx[N + 1];

            // Build z^16 - 1 (roots are 16th roots of unity)
            for (int i = 0; i <= N; i++) c[i] = 0.0;
            c[0] = -1.0;
            c[N] = 1.0;

            // Evaluate at a point near the unit circle
            var z = Cmplx.Polar(1.01, 0.37);

            PolySolver.PolynomialValDer(N, c, z, out _, out Cmplx derF, true);
            Cmplx derFD = FiniteDiffDer(N, c, z, true);
            AssertClose(derF, derFD, 1e-4, "degree 16 forward");

            var w = new Cmplx(0.5, 0.3);
            PolySolver.PolynomialValDer(N, c, w, out _, out Cmplx derB, false);
            Cmplx derFDB = FiniteDiffDer(N, c, w, false);
            AssertClose(derB, derFDB, 1e-4, "degree 16 backward");
        }

        // --- Ratio1 returns zero at an exact root ---

        [Fact]
        public void Ratio1_AtExactRoot_ReturnsZero()
        {
            // (z - 2)(z - 3) = z² - 5z + 6
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = 6.0;
            c[1] = -5.0;
            c[2] = 1.0;

            Cmplx r = PolySolver.Ratio1(2, c, 2.0, true);

            AssertClose(0.0, r, label: "ratio at root");
        }

        // --- Ratio1 gives a good Newton step toward a root ---

        [Fact]
        public void Ratio1_NearRoot_PointsTowardRoot()
        {
            // (z - 2)(z - 3) = z² - 5z + 6
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = 6.0;
            c[1] = -5.0;
            c[2] = 1.0;

            Cmplx z    = 2.1;
            Cmplx r    = PolySolver.Ratio1(2, c, z, true);
            Cmplx zNew = z - r;

            // Should be much closer to 2.0 than 2.1 was
            Assert.True((zNew - 2.0).Abs < (z - 2.0).Abs * 0.5,
                $"Newton step didn't improve: |z_new - 2| = {(zNew - 2.0).Abs}");
        }

        // --- PolynomialValDer2: second derivative checks ---

        [Fact]
        public void ValDer2_Constant_Forward()
        {
            Cmplx* c = stackalloc Cmplx[1];
            c[0] = new Cmplx(7.0, -3.0);

            PolySolver.PolynomialValDer2(0, c, new Cmplx(99.0),
                out Cmplx val, out Cmplx der, out Cmplx der2, true);

            AssertClose(c[0], val, label: "val");
            AssertClose(0.0, der, label: "der");
            AssertClose(0.0, der2, label: "der2");
        }

        [Fact]
        public void ValDer2_Linear_Forward()
        {
            // P(z) = 3 + 5z, P'(z) = 5, P''(z) = 0
            Cmplx* c = stackalloc Cmplx[2];
            c[0] = 3.0;
            c[1] = 5.0;

            PolySolver.PolynomialValDer2(1, c, 2.0,
                out Cmplx val, out Cmplx der, out Cmplx der2, true);

            AssertClose(13.0, val, label: "val");
            AssertClose(5.0, der, label: "der");
            AssertClose(0.0, der2, label: "der2");
        }

        [Fact]
        public void ValDer2_Quadratic_Forward()
        {
            // P(z) = 1 + 2z + 3z², P'(z) = 2 + 6z, P''(z) = 6
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = 1.0;
            c[1] = 2.0;
            c[2] = 3.0;
            Cmplx z = 4.0;

            PolySolver.PolynomialValDer2(2, c, z,
                out Cmplx val, out Cmplx der, out Cmplx der2, true);

            AssertClose(57.0, val, label: "val"); // 1 + 8 + 48
            AssertClose(26.0, der, label: "der"); // 2 + 24
            AssertClose(6.0, der2, label: "der2");
        }

        [Fact]
        public void ValDer2_Cubic_Forward()
        {
            // P(z) = 2 - z + 3z² + z³
            // P'(z) = -1 + 6z + 3z²
            // P''(z) = 6 + 6z
            Cmplx* c = stackalloc Cmplx[4];
            c[0] = 2.0;
            c[1] = -1.0;
            c[2] = 3.0;
            c[3] = 1.0;
            Cmplx z = 2.0;

            PolySolver.PolynomialValDer2(3, c, z,
                out Cmplx val, out Cmplx der, out Cmplx der2, true);

            AssertClose(20.0, val, label: "val");   // 2 - 2 + 12 + 8
            AssertClose(23.0, der, label: "der");   // -1 + 12 + 12
            AssertClose(18.0, der2, label: "der2"); // 6 + 12
        }

        [Fact]
        public void ValDer2_Quadratic_Backward()
        {
            // Q(w) = c[2] + c[1]w + c[0]w² = 3 + 2w + w²
            // Q'(w) = 2 + 2w
            // Q''(w) = 2
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = 1.0;
            c[1] = 2.0;
            c[2] = 3.0;
            Cmplx w = 4.0;

            PolySolver.PolynomialValDer2(2, c, w,
                out Cmplx val, out Cmplx der, out Cmplx der2, false);

            AssertClose(27.0, val, label: "val"); // 3 + 8 + 16
            AssertClose(10.0, der, label: "der"); // 2 + 8
            AssertClose(2.0, der2, label: "der2");
        }

        [Fact]
        public void ValDer2_MatchesValDer_ForFirstTwoOutputs()
        {
            // Position and first derivative should match the 2-output version
            Cmplx* c = stackalloc Cmplx[6];
            c[0] = new Cmplx(1.0, -2.0);
            c[1] = new Cmplx(3.0, 1.0);
            c[2] = new Cmplx(-1.0, 4.0);
            c[3] = new Cmplx(2.0, -3.0);
            c[4] = new Cmplx(0.5, 0.5);
            c[5] = new Cmplx(-1.0, 1.0);
            var z = new Cmplx(0.7, -0.3);

            PolySolver.PolynomialValDer(5, c, z, out Cmplx val1, out Cmplx der1, true);
            PolySolver.PolynomialValDer2(5, c, z, out Cmplx val2, out Cmplx der2, out _, true);

            AssertClose(val1, val2, label: "val forward");
            AssertClose(der1, der2, label: "der forward");

            var w = new Cmplx(0.4, 0.6);
            PolySolver.PolynomialValDer(5, c, w, out val1, out der1, false);
            PolySolver.PolynomialValDer2(5, c, w, out val2, out der2, out _, false);

            AssertClose(val1, val2, label: "val backward");
            AssertClose(der1, der2, label: "der backward");
        }

        [Fact]
        public void ValDer2_SecondDerivMatchesFiniteDifference_Forward()
        {
            Cmplx* c = stackalloc Cmplx[6];
            c[0] = new Cmplx(1.0, -2.0);
            c[1] = new Cmplx(3.0, 1.0);
            c[2] = new Cmplx(-1.0, 4.0);
            c[3] = new Cmplx(2.0, -3.0);
            c[4] = new Cmplx(0.5, 0.5);
            c[5] = new Cmplx(-1.0, 1.0);
            var z = new Cmplx(0.7, -0.3);

            PolySolver.PolynomialValDer2(5, c, z, out _, out _, out Cmplx der2, true);

            // FD of first derivative
            double h = 1e-7;
            PolySolver.PolynomialValDer(5, c, z + h, out _, out Cmplx dp, true);
            PolySolver.PolynomialValDer(5, c, z - h, out _, out Cmplx dm, true);
            Cmplx der2FD = (dp - dm) / (2.0 * h);

            AssertClose(der2, der2FD, 1e-4, "der2 forward vs FD");
        }

        [Fact]
        public void ValDer2_SecondDerivMatchesFiniteDifference_Backward()
        {
            Cmplx* c = stackalloc Cmplx[6];
            c[0] = new Cmplx(1.0, -2.0);
            c[1] = new Cmplx(3.0, 1.0);
            c[2] = new Cmplx(-1.0, 4.0);
            c[3] = new Cmplx(2.0, -3.0);
            c[4] = new Cmplx(0.5, 0.5);
            c[5] = new Cmplx(-1.0, 1.0);
            var w = new Cmplx(0.4, 0.6);

            PolySolver.PolynomialValDer2(5, c, w, out _, out _, out Cmplx der2, false);

            double h = 1e-7;
            PolySolver.PolynomialValDer(5, c, w + h, out _, out Cmplx dp, false);
            PolySolver.PolynomialValDer(5, c, w - h, out _, out Cmplx dm, false);
            Cmplx der2FD = (dp - dm) / (2.0 * h);

            AssertClose(der2, der2FD, 1e-4, "der2 backward vs FD");
        }

        [Fact]
        public void ValDer2_AtDoubleRoot_Der2IsNonzero()
        {
            // (z - 2)²(z + 1) = z³ - 3z² + 4
            // P''(2) = 6·2 - 6 = 6
            Cmplx* c = stackalloc Cmplx[4];
            c[0] = 4.0;
            c[1] = 0.0;
            c[2] = -3.0;
            c[3] = 1.0;

            PolySolver.PolynomialValDer2(3, c, 2.0,
                out Cmplx val, out Cmplx der, out Cmplx der2, true);

            AssertClose(0.0, val, label: "val at double root");
            AssertClose(0.0, der, label: "der at double root");
            AssertClose(6.0, der2, label: "der2 at double root");
        }

        [Fact]
        public void ValDer2_Degree16_FiniteDifferenceCheck()
        {
            const int N = 16;
            Cmplx*    c = stackalloc Cmplx[N + 1];

            // z^16 - 1
            for (int i = 0; i <= N; i++) c[i] = 0.0;
            c[0] = -1.0;
            c[N] = 1.0;

            var z = Cmplx.Polar(1.01, 0.37);
            PolySolver.PolynomialValDer2(N, c, z, out _, out _, out Cmplx der2, true);

            double h = 1e-7;
            PolySolver.PolynomialValDer(N, c, z + h, out _, out Cmplx dp, true);
            PolySolver.PolynomialValDer(N, c, z - h, out _, out Cmplx dm, true);
            Cmplx der2FD = (dp - dm) / (2.0 * h);

            AssertClose(der2, der2FD, 1e-2, "degree 16 der2 forward");
        }

        [Fact]
        public void ValDer2_ComplexCoefficients()
        {
            // P(z) = (1+i) + (2-i)z + (3+2i)z²
            // P''(z) = 2·(3+2i) = 6+4i (constant)
            Cmplx* c = stackalloc Cmplx[3];
            c[0] = new Cmplx(1.0, 1.0);
            c[1] = new Cmplx(2.0, -1.0);
            c[2] = new Cmplx(3.0, 2.0);

            PolySolver.PolynomialValDer2(2, c, new Cmplx(1.0, 1.0),
                out _, out _, out Cmplx der2, true);

            AssertClose(new Cmplx(6.0, 4.0), der2, label: "der2 complex coefficients");
        }
    }
}
