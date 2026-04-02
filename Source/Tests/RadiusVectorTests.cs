using System;
using Xunit;

namespace PatchedConicFixes.Tests
{
    public unsafe class RadiusVectorTests
    {
        const double Tol = 1e-10;

        static void AssertClose(double expected, double actual, double tol = Tol, string label = "")
        {
            Assert.True(Math.Abs(expected - actual) < tol,
                $"{label} Expected {expected:G17}, got {actual:G17}, diff = {Math.Abs(expected - actual):G17}");
        }

        static double Dot(double* a, double* b)
            => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];

        static double Norm(double* a)
            => Math.Sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);

        // --- RadiusVector: basic geometry ---

        [Fact]
        public void CircularEquatorial_Periapsis()
        {
            // Circular orbit (e=0), i=0, w=0, Om=0, u=0
            // At periapsis: r = a*(cos(0)-0)·P + a*sin(0)·η·Q = a·P = a·x̂
            var O = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, data.a1, data.e1, 0.0, r);

            AssertClose(2.0, r[0], label: "r[0]");
            AssertClose(0.0, r[1], label: "r[1]");
            AssertClose(0.0, r[2], label: "r[2]");
        }

        [Fact]
        public void CircularEquatorial_90Degrees()
        {
            // u=π/2: r = a*(cos(π/2)-0)·P + a*sin(π/2)·1·Q = a·Q = a·ŷ
            var O = new COrbitData(3.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, data.a1, data.e1, Math.PI / 2.0, r);

            AssertClose(0.0, r[0], label: "r[0]");
            AssertClose(3.0, r[1], label: "r[1]");
            AssertClose(0.0, r[2], label: "r[2]");
        }

        [Fact]
        public void CircularEquatorial_180Degrees()
        {
            // u=π: r = a*(cos(π)-0)·P = -a·P
            var O = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, data.a1, data.e1, Math.PI, r);

            AssertClose(-2.0, r[0], label: "r[0]");
            AssertClose(0.0, r[1], label: "r[1]");
            AssertClose(0.0, r[2], label: "r[2]");
        }

        [Fact]
        public void EllipticOrbit_RadiusMagnitude()
        {
            // |r| = a(1 - e·cos u)
            double a = 2.0, e = 0.3, u = 1.0;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, data.a1, data.e1, u, r);

            double expectedMag = a * (1.0 - e * Math.Cos(u));
            AssertClose(expectedMag, Norm(r), label: "|r|");
        }

        [Fact]
        public void EllipticOrbit_LiesInOrbitalPlane()
        {
            // r should be perpendicular to the orbital normal P × Q
            var O = new COrbitData(2.0, 0.3, 0.8, 1.2, 0.5);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, data.a1, data.e1, 1.7, r);

            // n = P × Q
            double nx = data.P1[1] * data.Q1[2] - data.P1[2] * data.Q1[1];
            double ny = data.P1[2] * data.Q1[0] - data.P1[0] * data.Q1[2];
            double nz = data.P1[0] * data.Q1[1] - data.P1[1] * data.Q1[0];

            double rdotn = r[0] * nx + r[1] * ny + r[2] * nz;
            AssertClose(0.0, rdotn, label: "r·n (in-plane check)");
        }

        [Fact]
        public void EllipticOrbit_Periapsis_AlongP()
        {
            // At u=0: r = a(1-e)·P (periapsis distance along P direction)
            double a = 3.0, e = 0.5;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, data.a1, data.e1, 0.0, r);

            double rp = a * (1.0 - e);
            // r should equal rp * P̂
            AssertClose(rp * data.P1[0], r[0], label: "r[0]");
            AssertClose(rp * data.P1[1], r[1], label: "r[1]");
            AssertClose(rp * data.P1[2], r[2], label: "r[2]");
        }

        [Fact]
        public void EllipticOrbit_Apoapsis_AlongNegP()
        {
            // At u=π: r = a(1+e)·(-P̂) = -a(1+e)·P̂
            double a = 3.0, e = 0.5;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, data.a1, data.e1, Math.PI, r);

            double ra = a * (1.0 + e);
            AssertClose(-ra * data.P1[0], r[0], label: "r[0]");
            AssertClose(-ra * data.P1[1], r[1], label: "r[1]");
            AssertClose(-ra * data.P1[2], r[2], label: "r[2]");
        }

        // --- RadiusVectorValDer: derivative checks ---

        [Fact]
        public void ValDer_PositionMatchesRadiusVector()
        {
            double a = 2.0, e = 0.3, u = 1.5;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r1 = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, a, e, u, r1);

            double* r2 = stackalloc double[3];
            double* rd = stackalloc double[3];
            Baluev.RadiusVectorValDer(true, data.P1, data.Q1, a, e, u, r2, rd);

            AssertClose(r1[0], r2[0], label: "r[0]");
            AssertClose(r1[1], r2[1], label: "r[1]");
            AssertClose(r1[2], r2[2], label: "r[2]");
        }

        [Fact]
        public void ValDer_DerivativeMatchesFiniteDifference()
        {
            double a = 2.0, e = 0.3, u = 1.5;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            double* rd = stackalloc double[3];
            Baluev.RadiusVectorValDer(true, data.P1, data.Q1, a, e, u, r, rd);

            // Finite difference: (r(u+h) - r(u-h)) / (2h) / a
            double h = 1e-7;
            double* rp = stackalloc double[3];
            double* rm = stackalloc double[3];
            Baluev.RadiusVector(true, data.P1, data.Q1, a, e, u + h, rp);
            Baluev.RadiusVector(true, data.P1, data.Q1, a, e, u - h, rm);

            for (int i = 0; i < 3; i++)
            {
                double fdDer = (rp[i] - rm[i]) / (2.0 * h * a);
                AssertClose(fdDer, rd[i], 1e-5, $"rd[{i}]");
            }
        }

        [Fact]
        public void ValDer_DerivativePerpToRadiusForCircular()
        {
            // For e=0: dr/du is perpendicular to r (uniform circular motion)
            double a = 2.0, e = 0.0, u = 1.0;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            double* rd = stackalloc double[3];
            Baluev.RadiusVectorValDer(true, data.P1, data.Q1, a, e, u, r, rd);

            // rd is dr/du / a, so actual derivative is a*rd
            // r · (a*rd) should be zero for circular orbit
            double dotRRd = r[0] * rd[0] + r[1] * rd[1] + r[2] * rd[2];
            AssertClose(0.0, dotRRd * a, 1e-10, "r · dr/du");
        }

        // --- RadiusVectorValDer2: second derivative checks ---

        [Fact]
        public void ValDer2_PositionAndFirstDerivMatchValDer()
        {
            double a = 2.0, e = 0.3, u = 1.5;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r1 = stackalloc double[3];
            double* rd1 = stackalloc double[3];
            Baluev.RadiusVectorValDer(true, data.P1, data.Q1, a, e, u, r1, rd1);

            double* r2 = stackalloc double[3];
            double* rd2 = stackalloc double[3];
            double* rdd = stackalloc double[3];
            Baluev.RadiusVectorValDer2(true, data.P1, data.Q1, a, e, u, r2, rd2, rdd);

            for (int i = 0; i < 3; i++)
            {
                AssertClose(r1[i], r2[i], label: $"r[{i}]");
                AssertClose(rd1[i], rd2[i], label: $"rd[{i}]");
            }
        }

        [Fact]
        public void ValDer2_SecondDerivMatchesFiniteDifference()
        {
            double a = 2.0, e = 0.3, u = 1.5;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            double* rd = stackalloc double[3];
            double* rdd = stackalloc double[3];
            Baluev.RadiusVectorValDer2(true, data.P1, data.Q1, a, e, u, r, rd, rdd);

            // FD of first derivative: (rd(u+h) - rd(u-h)) / (2h)
            double h = 1e-7;
            double* rp = stackalloc double[3];
            double* rdp = stackalloc double[3];
            double* rm = stackalloc double[3];
            double* rdm = stackalloc double[3];
            Baluev.RadiusVectorValDer(true, data.P1, data.Q1, a, e, u + h, rp, rdp);
            Baluev.RadiusVectorValDer(true, data.P1, data.Q1, a, e, u - h, rm, rdm);

            for (int i = 0; i < 3; i++)
            {
                double fdDer2 = (rdp[i] - rdm[i]) / (2.0 * h);
                AssertClose(fdDer2, rdd[i], 1e-5, $"rdd[{i}]");
            }
        }

        [Fact]
        public void ValDer2_CircularOrbit_SecondDerivProportionalToPosition()
        {
            // For circular e=0: d²r/du² = -r/a (since r = a·(cos u·P + sin u·Q))
            // rdd = d²r/du² / a, so rdd = -r/a²
            double a = 3.0, e = 0.0, u = 1.2;
            var O = new COrbitData(a, e, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double* r = stackalloc double[3];
            double* rd = stackalloc double[3];
            double* rdd = stackalloc double[3];
            Baluev.RadiusVectorValDer2(true, data.P1, data.Q1, a, e, u, r, rd, rdd);

            for (int i = 0; i < 3; i++)
                AssertClose(-r[i] / a, rdd[i], 1e-10, $"rdd[{i}] = -r[{i}]/a");
        }

        // --- DistanceBetween ---

        [Fact]
        public void DistanceBetween_SamePoint_IsZero()
        {
            var O = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var data = SAuxData.Create(in O, in O);

            double dist = Baluev.DistanceBetween(&data, true, true, 1.0, 1.0);
            AssertClose(0.0, dist, 1e-10, "same point distance");
        }

        [Fact]
        public void DistanceBetween_OppositePointsOnCircle()
        {
            // Two identical circular orbits, u1=0, u2=π → distance = 2a
            var O = new COrbitData(3.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O, in O);

            double dist = Baluev.DistanceBetween(&data, true, true, 0.0, Math.PI);
            AssertClose(6.0, dist, label: "diameter distance");
        }

        [Fact]
        public void DistanceBetween_90DegreesOnCircle()
        {
            // Identical circular orbit, u1=0, u2=π/2 → distance = a√2
            var O = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O, in O);

            double dist = Baluev.DistanceBetween(&data, true, true, 0.0, Math.PI / 2.0);
            AssertClose(2.0 * Math.Sqrt(2.0), dist, label: "90° distance");
        }

        [Fact]
        public void DistanceBetween_ConcentricCircles_SameAngle()
        {
            // Two coplanar circular orbits of different radii at the same angle
            // distance = |a1 - a2|
            var O1 = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var O2 = new COrbitData(5.0, 0.0, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            double dist = Baluev.DistanceBetween(&data, true, true, 0.7, 0.7);
            AssertClose(3.0, dist, label: "concentric same angle");
        }

        [Fact]
        public void DistanceBetween_PerpendicularPlanes()
        {
            // O1 equatorial, O2 polar, both circular, same a
            // At u1=0, u2=0: both at periapsis along their P vectors
            // O1: r1 = a·(1,0,0), O2: r2 = a·(1,0,0) — same point for w=0, Om=0
            var O1 = new COrbitData(2.0, 0.0, 0.0, 0.0, 0.0);
            var O2 = new COrbitData(2.0, 0.0, Math.PI / 2.0, 0.0, 0.0);
            var data = SAuxData.Create(in O1, in O2);

            // At u=0 both have r = a·P, and P for both starts at (1,0,0)
            double dist0 = Baluev.DistanceBetween(&data, true, true, 0.0, 0.0);
            AssertClose(0.0, dist0, 1e-10, "both at periapsis");

            // At u1=π/2, u2=π/2:
            // O1: r = a·(0,1,0), O2: r = a·(0,0,1)
            double dist90 = Baluev.DistanceBetween(&data, true, true, Math.PI / 2.0, Math.PI / 2.0);
            AssertClose(2.0 * Math.Sqrt(2.0), dist90, label: "90° on perpendicular planes");
        }

        [Fact]
        public void DistanceBetween_Symmetric()
        {
            // d(O1@u1, O2@u2) should equal d(O2@u2, O1@u1)
            var O1 = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data12 = SAuxData.Create(in O1, in O2);
            var data21 = SAuxData.Create(in O2, in O1);

            double d12 = Baluev.DistanceBetween(&data12, true, true, 1.0, 2.0);
            double d21 = Baluev.DistanceBetween(&data21, true, true, 2.0, 1.0);
            AssertClose(d12, d21, label: "distance symmetry");
        }

        [Fact]
        public void DistanceBetween_TriangleInequality()
        {
            var O1 = new COrbitData(2.0, 0.3, 0.5, 0.8, 1.2);
            var O2 = new COrbitData(3.0, 0.5, 0.7, 0.2, 1.5);
            var data = SAuxData.Create(in O1, in O2);

            // Three points: A = O1@0, B = O2@1, C = O2@2
            // d(A,C) ≤ d(A,B) + d(B,C)
            // For d(B,C) we need same-orbit data
            var data22 = SAuxData.Create(in O2, in O2);

            double dAB = Baluev.DistanceBetween(&data, true, true, 0.0, 1.0);
            double dAC = Baluev.DistanceBetween(&data, true, true, 0.0, 2.0);
            double dBC = Baluev.DistanceBetween(&data22, true, true, 1.0, 2.0);

            Assert.True(dAC <= dAB + dBC + Tol,
                $"Triangle inequality violated: {dAC:G17} > {dAB:G17} + {dBC:G17}");
        }

        [Fact]
        public void DistanceBetween_EllipticOrbit_KnownGeometry()
        {
            // Periapsis to apoapsis on same elliptic orbit:
            // d = a(1+e) - (-a(1-e)) along the apse line... wait, they're collinear
            // d = a(1+e) + a(1-e) = 2a
            double a = 3.0, e = 0.5;
            var O = new COrbitData(a, e, 0.0, 0.0, 0.0);
            var data = SAuxData.Create(in O, in O);

            double dist = Baluev.DistanceBetween(&data, true, true, 0.0, Math.PI);
            AssertClose(2.0 * a, dist, label: "periapsis to apoapsis");
        }
    }
}
