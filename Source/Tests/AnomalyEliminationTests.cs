using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public unsafe class AnomalyEliminationTests
    {
        const double Tol = 1e-10;

        static void AssertClose(double expected, double actual, double tol = Tol, string label = "")
        {
            Assert.True(Math.Abs(expected - actual) < tol,
                $"{label} Expected {expected:G17}, got {actual:G17}, diff = {Math.Abs(expected - actual):G17}");
        }

        // --- EliminatedAnomaly: basic geometry ---

        [Fact]
        public void IdenticalCircularOrbits_U2EqualsU1()
        {
            // Same circular orbit: the closest point at any u1 is u2 = u1
            var O = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O, in O);

            double u1 = 1.0;
            bool neg = Baluev.EliminatedAnomaly(&data, u1, out double u2, out double u2Alt);

            Assert.False(neg, "Should not have negative discriminant");
            // One of the candidates should be u1
            double err0 = Math.Abs(Baluev.AngleWrap(u2 - u1));
            double err1 = Math.Abs(Baluev.AngleWrap(u2Alt - u1));
            Assert.True(Math.Min(err0, err1) < 1e-8,
                $"Neither u2={u2:G17} nor u2Alt={u2Alt:G17} matches u1={u1:G17}");
        }

        [Fact]
        public void ConcentricCoplanarCircles_U2EqualsU1()
        {
            // Two coplanar circles, different radii: closest point is at same angle
            var O1 = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var O2 = new COrbitData(5.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            double u1 = 0.7;
            bool neg = Baluev.EliminatedAnomaly(&data, u1, out double u2, out double u2Alt);

            Assert.False(neg);
            double err0 = Math.Abs(Baluev.AngleWrap(u2 - u1));
            double err1 = Math.Abs(Baluev.AngleWrap(u2Alt - u1));
            Assert.True(Math.Min(err0, err1) < 1e-8,
                $"Neither u2={u2:G17} nor u2Alt={u2Alt:G17} matches u1={u1:G17}");
        }

        [Fact]
        public void EliminatedAnomaly_ProducesTwoCandidates()
        {
            // General configuration: two candidates should be distinct
            var O1 = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            bool neg = Baluev.EliminatedAnomaly(&data, 1.0, out double u2, out double u2Alt);

            Assert.False(neg);
            Assert.True(Math.Abs(Baluev.AngleWrap(u2 - u2Alt)) > 1e-6,
                $"Two candidates should be distinct: u2={u2:G17}, u2Alt={u2Alt:G17}");
        }

        [Fact]
        public void EliminatedAnomaly_EllipticHighEccentricity()
        {
            var O1 = new COrbitData(2.0, 0.9, 0.3, 0.5, 0.0);
            var O2 = new COrbitData(3.0, 0.1, 0.6, 1.0, 0.8);
            var data = SAuxData.Create(in O1, in O2);

            double u1 = 0.5;
            bool neg = Baluev.EliminatedAnomaly(&data, u1, out double u2, out _);

            // Even if neg, the returned u2 should give a finite distance
            double d = Baluev.DistanceBetween(&data, true, true, u1, u2);
            Assert.True(!double.IsNaN(d) && !double.IsInfinity(d), $"Distance should be finite, got {d}");
            Assert.True(d >= 0, $"Distance should be non-negative, got {d}");
        }

        // --- MinDistanceFor ---

        [Fact]
        public void MinDistanceFor_ConcentricCircles()
        {
            // Two coplanar circles: min distance = |a1 - a2| at u2 = u1
            var O1 = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var O2 = new COrbitData(5.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            double u1 = 0.7;
            double d = Baluev.MinDistanceFor(&data, u1, out double u2);

            AssertClose(3.0, d, 1e-8, "distance");
            AssertClose(u1, u2, 1e-8, "u2 = u1");
        }

        [Fact]
        public void MinDistanceFor_IdenticalCircularOrbits()
        {
            var O = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O, in O);

            double d = Baluev.MinDistanceFor(&data, 1.0, out double u2);

            AssertClose(0.0, Math.Abs(d), 1e-8, "distance on same orbit");
        }

        [Fact]
        public void MinDistanceFor_ConsistentWithDistanceBetween()
        {
            // The returned distance should match DistanceBetween(u1, u2)
            var O1 = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            double u1 = 1.0;
            double d = Baluev.MinDistanceFor(&data, u1, out double u2);

            double dCheck = Baluev.DistanceBetween(&data, true, true, u1, u2);
            AssertClose(Math.Abs(d), dCheck, 1e-10, "distance consistency");
        }

        [Fact]
        public void MinDistanceFor_ReturnsSmaller()
        {
            // The returned distance should be ≤ distance at the alternative u2
            var O1 = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            double u1 = 1.0;
            Baluev.EliminatedAnomaly(&data, u1, out double u2a, out double u2b);
            double da = Baluev.DistanceBetween(&data, true, true, u1, u2a);
            double db = Baluev.DistanceBetween(&data, true, true, u1, u2b);

            double d = Baluev.MinDistanceFor(&data, u1, out _);
            double expectedMin = Math.Min(da, db);
            AssertClose(expectedMin, Math.Abs(d), 1e-10, "returns min");
        }

        [Fact]
        public void MinDistanceFor_NegativeReturnMeansDegenerate()
        {
            // Find a case where the discriminant is negative
            // Identical orbit at periapsis: this is actually well-behaved
            // Use near-coplanar orbits where the geometry degenerates
            var O1 = new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0);
            var O2 = new COrbitData(1.0, 0.0, 1e-10, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            // Sweep for a degenerate case
            for (int k = 0; k < 32; k++)
            {
                double u1 = -Math.PI + k * 2.0 * Math.PI / 32.0;
                double d = Baluev.MinDistanceFor(&data, u1, out _);
                if (d < 0)
                {
                    // Absolute value should still be a valid distance
                    Assert.True(!double.IsNaN(Math.Abs(d)));
                    Assert.True(!double.IsInfinity(Math.Abs(d)));
                    break;
                }
            }
            // Note: not all configurations produce negative discriminant;
            // this test just verifies the code path doesn't crash
        }

        [Fact]
        public void MinDistanceFor_SeveralU1Values()
        {
            // Sweep u1, verify consistency at each point
            var O1 = new COrbitData(1.5, 0.2, 0.3, 0.5, 0.0);
            var O2 = new COrbitData(2.5, 0.4, 0.6, 1.0, 0.8);
            var data = SAuxData.Create(in O1, in O2);

            for (int k = 0; k < 16; k++)
            {
                double u1 = -Math.PI + k * 2.0 * Math.PI / 16.0;
                double d = Baluev.MinDistanceFor(&data, u1, out double u2);

                if (d >= 0)
                {
                    double dCheck = Baluev.DistanceBetween(&data, true, true, u1, u2);
                    AssertClose(d, dCheck, 1e-8, $"consistency at u1={u1:F3}");
                }
                else
                {
                    // Negative return: absolute value should still be close
                    double dCheck = Baluev.DistanceBetween(&data, true, true, u1, u2);
                    AssertClose(-d, dCheck, 1e-8, $"neg consistency at u1={u1:F3}");
                }
            }
        }

        [Fact]
        public void MinDistanceFor_SymmetricOrbits()
        {
            // Two orbits with same shape but different w:
            // Sweeping u1 should give symmetric distance profile
            var O1 = new COrbitData(2.0, 0.3, 0.5, 0.0, 0.0);
            var O2 = new COrbitData(2.0, 0.3, 0.5, Math.PI, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            // At u1=0 and u1=π, the geometry should be related by symmetry
            double d0 = Math.Abs(Baluev.MinDistanceFor(&data, 0.0, out _));
            double dPi = Math.Abs(Baluev.MinDistanceFor(&data, Math.PI, out _));
            AssertClose(d0, dPi, 1e-8, "symmetric distance");
        }
    }
}
