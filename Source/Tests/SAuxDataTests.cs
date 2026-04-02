using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public class SAuxDataTests
    {
        private const double Tol = 1e-12;

        private static void AssertClose(double expected, double actual, double tol = Tol, string label = "")
        {
            Assert.True(Math.Abs(expected - actual) < tol,
                $"{label} Expected {expected:G17}, got {actual:G17}, diff = {Math.Abs(expected - actual):G17}");
        }

        // --- COrbitData P/Q vector sanity ---

        [Fact]
        public void COrbitData_PQ_AreOrthonormal()
        {
            // Arbitrary orbit
            var O = new COrbitData(2.0, 0.3, 0.5, 1.2, 0.8);

            double PdotP = O.P0 * O.P0 + O.P1 * O.P1 + O.P2 * O.P2;
            double QdotQ = O.Q0 * O.Q0 + O.Q1 * O.Q1 + O.Q2 * O.Q2;
            double PdotQ = O.P0 * O.Q0 + O.P1 * O.Q1 + O.P2 * O.Q2;

            AssertClose(1.0, PdotP, label: "|P|²");
            AssertClose(1.0, QdotQ, label: "|Q|²");
            AssertClose(0.0, PdotQ, label: "P·Q");
        }

        [Fact]
        public void COrbitData_ZeroAngles_PIsXhat_QIsYhat()
        {
            // i=0, w=0, Om=0 → P = x̂, Q = ŷ
            var O = new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0);

            AssertClose(1.0, O.P0, label: "P[0]");
            AssertClose(0.0, O.P1, label: "P[1]");
            AssertClose(0.0, O.P2, label: "P[2]");
            AssertClose(0.0, O.Q0, label: "Q[0]");
            AssertClose(1.0, O.Q1, label: "Q[1]");
            AssertClose(0.0, O.Q2, label: "Q[2]");
        }

        [Fact]
        public void COrbitData_W90_RotatesP()
        {
            // i=0, Om=0, w=π/2 → P should point along ŷ, Q along -x̂
            var O = new COrbitData(1.0, 0.0, 0.0, Math.PI / 2.0, 0.0);

            AssertClose(0.0, O.P0, 1e-10, "P[0]");
            AssertClose(1.0, O.P1, 1e-10, "P[1]");
            AssertClose(0.0, O.P2, 1e-10, "P[2]");

            AssertClose(-1.0, O.Q0, 1e-10, "Q[0]");
            AssertClose(0.0, O.Q1, 1e-10, "Q[1]");
            AssertClose(0.0, O.Q2, 1e-10, "Q[2]");

            // Still orthonormal
            double PdotQ = O.P0 * O.Q0 + O.P1 * O.Q1 + O.P2 * O.Q2;
            AssertClose(0.0, PdotQ, 1e-10, "P·Q");
        }

        // --- identical orbits ---

        [Fact]
        public void IdenticalOrbits_DotProducts()
        {
            var O    = new COrbitData(2.0, 0.3, 0.5, 1.2, 0.8);
            var data = SAuxData.Create(in O, in O);

            // P1·P2 = |P|² = 1
            AssertClose(1.0, data.Pp, label: "Pp");

            // S·S' = (Q·η)·(Q·η) = η² · |Q|² = η² = 1-e²
            double eta = Math.Sqrt(1 - 0.3 * 0.3);
            AssertClose(eta * eta, data.Ss, label: "Ss");

            // P·S' = P·(Q·η) = 0 (P⊥Q)
            AssertClose(0.0, data.Ps, label: "Ps");
            AssertClose(0.0, data.Sp, label: "Sp");
        }

        [Fact]
        public void IdenticalOrbits_MutualInclination()
        {
            var O    = new COrbitData(2.0, 0.3, 0.5, 1.2, 0.8);
            var data = SAuxData.Create(in O, in O);

            AssertClose(0.0, data.I, 1e-10, "I");
        }

        [Fact]
        public void IdenticalOrbits_AlphaAndK()
        {
            var O    = new COrbitData(2.0, 0.3, 0.5, 1.2, 0.8);
            var data = SAuxData.Create(in O, in O);

            AssertClose(1.0, data.alpha1, label: "alpha1");
            AssertClose(1.0, data.alpha2, label: "alpha2");
            AssertClose(0.3 * 0.3, data.K, label: "K");
        }

        // --- coplanar orbits (same plane, different shape) ---

        [Fact]
        public void CoplanarOrbits_ZeroMutualInclination()
        {
            // Same i and Om, different a, e, w
            var O1   = new COrbitData(1.0, 0.2, 0.5, 0.3, 0.8);
            var O2   = new COrbitData(2.0, 0.4, 0.5, 1.0, 0.8);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(0.0, data.I, 1e-10, "I");
            AssertClose(0.0, data.abs_w, 1e-10, "abs_w");
        }

        [Fact]
        public void CoplanarOrbits_DotProductRelations()
        {
            // Coplanar: P1, Q1, P2, Q2 all lie in the same plane
            // So Pp² + Sp² = 1 (unit vector P2 decomposed into P1, Q1 basis)
            // and Ps² + Ss² = η2² (unit vector Q2·η2 decomposed)
            var O1   = new COrbitData(1.0, 0.2, 0.5, 0.3, 0.8);
            var O2   = new COrbitData(2.0, 0.4, 0.5, 1.0, 0.8);
            var data = SAuxData.Create(in O1, in O2);

            double eta1 = Math.Sqrt(1 - 0.2 * 0.2);
            double eta2 = Math.Sqrt(1 - 0.4 * 0.4);

            // P2 decomposed into P1/Q1 frame: Pp² + (Sp/η1)² = 1
            AssertClose(1.0, data.Pp * data.Pp + data.Sp / eta1 * (data.Sp / eta1),
                1e-10, "Pp²+(Sp/η1)²");

            // Q2·η2 decomposed: (Ps/η2)²·η2² + (Ss/(η1·η2))²·η1²·η2² ... simpler:
            // Ps² + Ss² should equal η2² times (dot(P1,Q2)² + dot(Q1,Q2)²·η1²)
            // Actually just check: Ps²/η2² + Ss²/(η1²·η2²) = 1
            AssertClose(1.0,
                data.Ps * data.Ps / (eta2 * eta2) +
                data.Ss * data.Ss / (eta1 * eta1 * eta2 * eta2),
                1e-10, "Q2 decomposition");
        }

        // --- perpendicular orbital planes ---

        [Fact]
        public void PerpendicularPlanes_MutualInclination90()
        {
            // O1 in equatorial plane (i=0), O2 in polar plane (i=π/2)
            var O1   = new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0);
            var O2   = new COrbitData(1.0, 0.0, Math.PI / 2.0, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(Math.PI / 2.0, data.I, 1e-10, "I");
        }

        // --- opposite orbital planes ---

        [Fact]
        public void RetrogradeOrbit_MutualInclination180()
        {
            var O1   = new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0);
            var O2   = new COrbitData(1.0, 0.0, Math.PI, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(Math.PI, data.I, 1e-10, "I");
        }

        // --- semi-latus rectum ---

        [Fact]
        public void SemiLatusRectum()
        {
            var O1   = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(2.0 * (1 - 0.3 * 0.3), data.p1, label: "p1");
            AssertClose(3.0 * (1 - 0.5 * 0.5), data.p2, label: "p2");
        }

        // --- alpha ratios ---

        [Fact]
        public void AlphaRatios()
        {
            var O1   = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(2.0 / 3.0, data.alpha1, label: "alpha1");
            AssertClose(3.0 / 2.0, data.alpha2, label: "alpha2");
            AssertClose(1.0, data.alpha1 * data.alpha2, label: "α1·α2");
        }

        // --- K value ---

        [Fact]
        public void KValue()
        {
            var O1   = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(data.alpha2 * data.e2 * data.e2, data.K, label: "K");
        }

        // --- circular orbits: e=0 simplifies everything ---

        [Fact]
        public void CircularOrbits_EtaIsOne()
        {
            var O1   = new COrbitData(1.0, 0.0, 0.3, 0.0, 0.0);
            var O2   = new COrbitData(2.0, 0.0, 0.3, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            // η = √(1-e²) = 1, so Ss = Q1·Q2, no scaling
            double Q1Q2 = O1.Q0 * O2.Q0 + O1.Q1 * O2.Q1 + O1.Q2 * O2.Q2;
            AssertClose(Q1Q2, data.Ss, label: "Ss");

            AssertClose(0.0, data.K, label: "K (e=0)");
        }

        // --- negative eccentricity handling ---

        [Fact]
        public void NegativeEccentricity_FlipsToPositive()
        {
            // Negative e should be treated as positive with ω shifted by π
            var O1   = new COrbitData(2.0, -0.3, 0.5, 1.0, 0.8);
            var O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            Assert.True(data.e1 > 0, "e1 should be positive");
            AssertClose(0.3, data.e1, label: "e1");
        }

        [Fact]
        public void NegativeEccentricity_MatchesShiftedOmega()
        {
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);

            var O1neg   = new COrbitData(2.0, -0.3, 0.5, 1.0, 0.8);
            var dataNeg = SAuxData.Create(in O1neg, in O2);

            // e gets flipped positive
            AssertClose(0.3, dataNeg.e1, label: "e1 flipped");

            // w gets shifted by π
            AssertClose(1.0 + Math.PI, dataNeg.w1, label: "w1 shifted");

            // Derived scalars are consistent with the flipped e
            AssertClose(dataNeg.a1 * (1 - dataNeg.e1 * dataNeg.e1), dataNeg.p1, label: "p1");
            AssertClose(dataNeg.alpha2 * dataNeg.e2 * dataNeg.e2, dataNeg.K, label: "K");
        }

        // --- hyperbolic orbit (e > 1) ---

        [Fact]
        public void HyperbolicOrbit_NegativeSemiMajorAxis()
        {
            var O1   = new COrbitData(2.0, 1.5, 0.3, 0.5, 0.0);
            var O2   = new COrbitData(3.0, 0.5, 0.5, 1.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            // e1 > 1 ⟹ a1 should be negative
            Assert.True(data.a1 < 0, "a1 should be negative for hyperbolic");
            AssertClose(-2.0, data.a1, label: "a1");

            // p = a(1-e²) should still work: p = -2·(1-2.25) = -2·(-1.25) = 2.5
            AssertClose(2.0 * (1.5 * 1.5 - 1), data.p1, label: "p1 hyperbolic");
        }

        [Fact]
        public void HyperbolicOrbit_EtaUsesAbsValue()
        {
            // η = √(|1-e²|) for hyperbolic
            var O1   = new COrbitData(2.0, 1.5, 0.0, 0.0, 0.0);
            var O2   = new COrbitData(2.0, 1.5, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            // Same orbit: Ss = η²·|Q|² = |1-e²| = |1-2.25| = 1.25
            double eta2 = Math.Abs(1 - 1.5 * 1.5);
            AssertClose(eta2, data.Ss, 1e-10, "Ss hyperbolic");
        }

        // --- symmetry: swapping O1/O2 ---

        [Fact]
        public void SwappedOrbits_MutualInclinationUnchanged()
        {
            var O1 = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);

            var data12 = SAuxData.Create(in O1, in O2);
            var data21 = SAuxData.Create(in O2, in O1);

            AssertClose(data12.I, data21.I, 1e-10, "I symmetric");
        }

        [Fact]
        public void SwappedOrbits_DotProductsTranspose()
        {
            var O1 = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);

            var data12 = SAuxData.Create(in O1, in O2);
            var data21 = SAuxData.Create(in O2, in O1);

            // Pp is symmetric: P1·P2 = P2·P1
            AssertClose(data12.Pp, data21.Pp, 1e-10, "Pp symmetric");

            // Ps and Sp swap: P1·S2 becomes P2·S1
            AssertClose(data12.Ps, data21.Sp, 1e-10, "Ps/Sp swap");
            AssertClose(data12.Sp, data21.Ps, 1e-10, "Sp/Ps swap");

            // Ss is symmetric: S1·S2 = S2·S1
            AssertClose(data12.Ss, data21.Ss, 1e-10, "Ss symmetric");
        }

        [Fact]
        public void SwappedOrbits_AlphasInvert()
        {
            var O1 = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);

            var data12 = SAuxData.Create(in O1, in O2);
            var data21 = SAuxData.Create(in O2, in O1);

            AssertClose(data12.alpha1, data21.alpha2, 1e-10, "alpha1/alpha2 swap");
            AssertClose(data12.alpha2, data21.alpha1, 1e-10, "alpha2/alpha1 swap");
        }

        // --- dot product cross-check via direct computation ---

        [Fact]
        public void DotProducts_CrossCheckDirect()
        {
            var O1   = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            double eta1 = Math.Sqrt(1 - 0.3 * 0.3);
            double eta2 = Math.Sqrt(1 - 0.5 * 0.5);

            double Pp   = O1.P0 * O2.P0 + O1.P1 * O2.P1 + O1.P2 * O2.P2;
            double PQ2  = O1.P0 * O2.Q0 + O1.P1 * O2.Q1 + O1.P2 * O2.Q2;
            double Q1P  = O1.Q0 * O2.P0 + O1.Q1 * O2.P1 + O1.Q2 * O2.P2;
            double Q1Q2 = O1.Q0 * O2.Q0 + O1.Q1 * O2.Q1 + O1.Q2 * O2.Q2;

            AssertClose(Pp, data.Pp, label: "Pp direct");
            AssertClose(PQ2 * eta2, data.Ps, label: "Ps direct");
            AssertClose(Q1P * eta1, data.Sp, label: "Sp direct");
            AssertClose(Q1Q2 * eta1 * eta2, data.Ss, label: "Ss direct");
        }

        // --- mutual inclination for known geometries ---

        [Fact]
        public void MutualInclination_45Degrees()
        {
            // O1 at i=0, O2 at i=π/4, same Om → I = π/4
            var O1   = new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0);
            var O2   = new COrbitData(1.0, 0.0, Math.PI / 4.0, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(Math.PI / 4.0, data.I, 1e-10, "I = 45°");
        }

        [Fact]
        public void MutualInclination_DifferentOm()
        {
            // Two equatorial orbits with different Om → same plane, I = 0
            var O1   = new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0);
            var O2   = new COrbitData(1.0, 0.0, 0.0, 0.0, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(0.0, data.I, 1e-10, "I equatorial");
        }

        [Fact]
        public void MutualInclination_ObtuseCosI()
        {
            // i1=0, i2=120° → cosI = cos(120°) < 0, I should be > π/2
            var O1   = new COrbitData(1.0, 0.0, 0.0, 0.0, 0.0);
            var O2   = new COrbitData(1.0, 0.0, 2.0 * Math.PI / 3.0, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(2.0 * Math.PI / 3.0, data.I, 1e-10, "I = 120°");
        }

        // --- w projections onto mutual node vector ---

        [Fact]
        public void NodeProjections_CoplanarOrbitsAreZero()
        {
            // Same plane ⟹ w = 0 ⟹ all projections zero
            var O1   = new COrbitData(1.0, 0.3, 0.5, 0.3, 0.8);
            var O2   = new COrbitData(2.0, 0.5, 0.5, 1.0, 0.8);
            var data = SAuxData.Create(in O1, in O2);

            AssertClose(0.0, data.P1w, 1e-10, "P1w");
            AssertClose(0.0, data.P2w, 1e-10, "P2w");
            AssertClose(0.0, data.Q1w, 1e-10, "Q1w");
            AssertClose(0.0, data.Q2w, 1e-10, "Q2w");
        }

        [Fact]
        public void NodeProjections_SumOfSquares()
        {
            // For each orbit: P_w² + Q_w² = |w_perp|² (projection of w into orbital plane)
            // This must be ≤ |w|²
            var O1   = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            double proj1 = data.P1w * data.P1w + data.Q1w * data.Q1w;
            double proj2 = data.P2w * data.P2w + data.Q2w * data.Q2w;
            double w2    = data.abs_w * data.abs_w;

            Assert.True(proj1 <= w2 + Tol,
                $"P1w²+Q1w² = {proj1:G17} exceeds |w|² = {w2:G17}");
            Assert.True(proj2 <= w2 + Tol,
                $"P2w²+Q2w² = {proj2:G17} exceeds |w|² = {w2:G17}");
        }

        // --- consistency: Pp² + Ps²/η2² + (out-of-plane)² = 1 ---

        [Fact]
        public void DotProductNormalization()
        {
            // P1 is a unit vector. Decompose into P2, Q2, and the normal to plane 2.
            // P1·P2 = Pp, P1·Q2 = Ps/η2, and the out-of-plane component completes the basis.
            // So Pp² + (Ps/η2)² + (P1·n2)² = 1
            var O1   = new COrbitData(2.0, 0.3, 0.5, 1.0, 0.8);
            var O2   = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            double eta2 = Math.Sqrt(1 - 0.5 * 0.5);

            // Normal to orbit 2: n2 = P2 × Q2
            double n2x = O2.P1 * O2.Q2 - O2.P2 * O2.Q1;
            double n2y = O2.P2 * O2.Q0 - O2.P0 * O2.Q2;
            double n2z = O2.P0 * O2.Q1 - O2.P1 * O2.Q0;

            double P1n2 = O1.P0 * n2x + O1.P1 * n2y + O1.P2 * n2z;

            double sum = data.Pp * data.Pp +
                data.Ps / eta2 * (data.Ps / eta2) +
                P1n2 * P1n2;

            AssertClose(1.0, sum, 1e-10, "P1 decomposition into orbit-2 frame");
        }
    }
}
