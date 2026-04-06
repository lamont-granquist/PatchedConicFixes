using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public class AngleHelperTests
    {
        private const double Tol = 1e-12;

        private static void AssertClose(double expected, double actual, double tol = Tol, string label = "")
        {
            Assert.True(Math.Abs(expected - actual) < tol,
                $"{label} Expected {expected:G17}, got {actual:G17}, diff = {Math.Abs(expected - actual):G17}");
        }

        // --- AngleWrap ---

        [Fact]
        public void AngleWrap_Zero() => AssertClose(0.0, Baluev.AngleWrap(0.0));

        [Fact]
        public void AngleWrap_AlreadyInRange()
        {
            AssertClose(1.0, Baluev.AngleWrap(1.0));
            AssertClose(-1.0, Baluev.AngleWrap(-1.0));
            AssertClose(0.5, Baluev.AngleWrap(0.5));
        }

        [Fact]
        public void AngleWrap_ExactlyPi()
        {
            // π should map to π or -π (both are valid boundary values)
            double result = Baluev.AngleWrap(Math.PI);
            Assert.True(Math.Abs(Math.Abs(result) - Math.PI) < Tol,
                $"AngleWrap(π) = {result:G17}, expected ±π");
        }

        [Fact]
        public void AngleWrap_ExactlyNegativePi()
        {
            double result = Baluev.AngleWrap(-Math.PI);
            Assert.True(Math.Abs(Math.Abs(result) - Math.PI) < Tol,
                $"AngleWrap(-π) = {result:G17}, expected ±π");
        }

        [Fact]
        public void AngleWrap_JustOverPi()
        {
            // π + 0.1 should wrap to roughly -π + 0.1
            double input    = Math.PI + 0.1;
            double expected = -Math.PI + 0.1;
            AssertClose(expected, Baluev.AngleWrap(input));
        }

        [Fact]
        public void AngleWrap_JustUnderNegativePi()
        {
            // -π - 0.1 should wrap to roughly π - 0.1
            double input    = -Math.PI - 0.1;
            double expected = Math.PI - 0.1;
            AssertClose(expected, Baluev.AngleWrap(input));
        }

        [Fact]
        public void AngleWrap_FullCircle()
        {
            AssertClose(0.0, Baluev.AngleWrap(2.0 * Math.PI), 1e-10);
            AssertClose(0.0, Baluev.AngleWrap(-2.0 * Math.PI), 1e-10);
        }

        [Fact]
        public void AngleWrap_MultipleRevolutions()
        {
            // 5π = 2·2π + π → should wrap to ±π
            double result = Baluev.AngleWrap(5.0 * Math.PI);
            Assert.True(Math.Abs(Math.Abs(result) - Math.PI) < 1e-10,
                $"AngleWrap(5π) = {result:G17}, expected ±π");

            // 4.5π = 2·2π + 0.5π → should wrap to 0.5π
            AssertClose(0.5 * Math.PI, Baluev.AngleWrap(4.5 * Math.PI), 1e-10);
        }

        [Fact]
        public void AngleWrap_NegativeMultipleRevolutions()
        {
            // -4.5π → should wrap to -0.5π
            AssertClose(-0.5 * Math.PI, Baluev.AngleWrap(-4.5 * Math.PI), 1e-10);
        }

        [Fact]
        public void AngleWrap_LargePositive()
        {
            // 100.3 radians — verify result is in [-π, π]
            double result = Baluev.AngleWrap(100.3);
            Assert.True(result >= -Math.PI && result <= Math.PI,
                $"AngleWrap(100.3) = {result:G17} out of [-π, π]");

            // Cross-check: sin and cos should match
            AssertClose(Math.Sin(100.3), Math.Sin(result), 1e-10, "sin");
            AssertClose(Math.Cos(100.3), Math.Cos(result), 1e-10, "cos");
        }

        [Fact]
        public void AngleWrap_LargeNegative()
        {
            double result = Baluev.AngleWrap(-100.3);
            Assert.True(result >= -Math.PI && result <= Math.PI,
                $"AngleWrap(-100.3) = {result:G17} out of [-π, π]");

            AssertClose(Math.Sin(-100.3), Math.Sin(result), 1e-10, "sin");
            AssertClose(Math.Cos(-100.3), Math.Cos(result), 1e-10, "cos");
        }

        // --- Atan2h ---

        [Fact]
        public void Atan2h_Zero()
        {
            // atanh(0/x) = 0 for any x ≠ 0
            AssertClose(0.0, Baluev.Atan2h(0.0, 1.0));
            AssertClose(0.0, Baluev.Atan2h(0.0, 5.0));
        }

        [Fact]
        public void Atan2h_MatchesAtanh()
        {
            // For x = 1: atan2h(y, 1) = atanh(y)
            // atanh(0.5) = (log(1.5) - log(0.5)) / 2
            double y        = 0.5;
            double expected = 0.5 * Math.Log(1.5 / 0.5); // atanh(0.5)
            AssertClose(expected, Baluev.Atan2h(y, 1.0));
        }

        [Fact]
        public void Atan2h_ScaleInvariance()
        {
            // atan2h(ky, kx) = atan2h(y, x) for k > 0
            double y    = 0.3, x = 2.0;
            double val1 = Baluev.Atan2h(y, x);
            double val2 = Baluev.Atan2h(3.0 * y, 3.0 * x);
            AssertClose(val1, val2, 1e-10, "scale invariance");
        }

        [Fact]
        public void Atan2h_OddInY()
        {
            // atan2h(-y, x) = -atan2h(y, x)
            double y = 0.7, x = 3.0;
            AssertClose(-Baluev.Atan2h(y, x), Baluev.Atan2h(-y, x), 1e-10, "odd in y");
        }

        [Fact]
        public void Atan2h_KnownValues()
        {
            // atanh(1/3) = atan2h(1, 3) = (log(4) - log(2)) / 2 = log(2)/2
            AssertClose(Math.Log(2.0) / 2.0, Baluev.Atan2h(1.0, 3.0));

            // atanh(1/2) = atan2h(1, 2) = (log(3) - log(1)) / 2 = log(3)/2
            AssertClose(Math.Log(3.0) / 2.0, Baluev.Atan2h(1.0, 2.0));
        }

        [Fact]
        public void Atan2h_SmallArgument()
        {
            // For small y/x, atanh(y/x) ≈ y/x
            double y = 1e-8, x = 1.0;
            AssertClose(y, Baluev.Atan2h(y, x), 1e-14, "small argument");
        }

        // --- Atan2Smart ---

        [Fact]
        public void Atan2Smart_Trig_MatchesAtan2()
        {
            double y = 3.0, x = 4.0;
            AssertClose(Math.Atan2(y, x), Baluev.Atan2Smart(y, x, true));
        }

        [Fact]
        public void Atan2Smart_Hyp_MatchesAtan2h()
        {
            double y = 0.3, x = 2.0;
            AssertClose(Baluev.Atan2h(y, x), Baluev.Atan2Smart(y, x, false));
        }

        [Fact]
        public void Atan2Smart_Trig_AllQuadrants()
        {
            // Just verify it delegates correctly across all quadrants
            AssertClose(Math.Atan2(1.0, 1.0), Baluev.Atan2Smart(1.0, 1.0, true), label: "Q1");
            AssertClose(Math.Atan2(1.0, -1.0), Baluev.Atan2Smart(1.0, -1.0, true), label: "Q2");
            AssertClose(Math.Atan2(-1.0, -1.0), Baluev.Atan2Smart(-1.0, -1.0, true), label: "Q3");
            AssertClose(Math.Atan2(-1.0, 1.0), Baluev.Atan2Smart(-1.0, 1.0, true), label: "Q4");
        }

        [Fact]
        public void Atan2Smart_Hyp_SignConsistency()
        {
            // Positive y → positive result
            Assert.True(Baluev.Atan2Smart(0.5, 2.0, false) > 0);
            // Negative y → negative result
            Assert.True(Baluev.Atan2Smart(-0.5, 2.0, false) < 0);
        }
    }
}
