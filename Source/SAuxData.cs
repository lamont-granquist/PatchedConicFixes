using Unity.Burst;
using Unity.Mathematics;

namespace PatchedConicFixes
{
    [BurstCompile]
    public unsafe struct SAuxData
    {
        // MOID computation
        public double e1,     e2, a1, a2;
        public double alpha1, alpha2;
        public double K;
        public double Pp, Ps, Sp, Ss;

        // Linking coefficients
        public double p1,  p2,  w1,  w2, I, abs_w;
        public double P1w, P2w, Q1w, Q2w;

        // Orbital basis vectors (inlined, no external pointers)
        public fixed double P1[3];
        public fixed double P2[3];
        public fixed double Q1[3];
        public fixed double Q2[3];

        public static SAuxData Create(in COrbitData O1, in COrbitData O2)
        {
            var d = new SAuxData();

            d.e1 = O1.e;
            d.e2 = O2.e;
            d.a1 = O1.a;
            d.a2 = O2.a;
            d.w1 = O1.w;
            d.w2 = O2.w;

            if (d.e1 < 0)
            {
                d.e1 =  -d.e1;
                d.w1 += math.PI_DBL;
            }

            if (d.e2 < 0)
            {
                d.e2 =  -d.e2;
                d.w2 += math.PI_DBL;
            }

            d.a1 = math.abs(d.a1);
            d.a2 = math.abs(d.a2);
            if (!(d.e1 <= 1)) d.a1 = -d.a1;
            if (!(d.e2 <= 1)) d.a2 = -d.a2;

            d.p1     = d.a1 * (1 - d.e1 * d.e1);
            d.p2     = d.a2 * (1 - d.e2 * d.e2);
            d.alpha1 = d.a1 / d.a2;
            d.alpha2 = d.a2 / d.a1;
            d.K      = d.alpha2 * d.e2 * d.e2;

            double eta1 = math.sqrt(math.abs(1 - d.e1 * d.e1));
            double eta2 = math.sqrt(math.abs(1 - d.e2 * d.e2));

            double i1  = O1.i,  i2  = O2.i;
            double Om1 = O1.Om, Om2 = O2.Om;

            double c1 = math.cos(i1), s1 = math.sin(i1);
            double c2 = math.cos(i2), s2 = math.sin(i2);

            // Mutual node vector
            double w0 = c1 * s2 * math.cos(Om2) - s1 * c2 * math.cos(Om1);
            double w1 = c1 * s2 * math.sin(Om2) - s1 * c2 * math.sin(Om1);
            double w2 = s1 * s2 * math.sin(Om2 - Om1);

            d.abs_w = math.sqrt(w0 * w0 + w1 * w1 + w2 * w2);

            double cosI         = c1 * c2 + s1 * s2 * math.cos(Om2 - Om1);
            double absW_clamped = math.clamp(d.abs_w, -1.0, 1.0);
            d.I = cosI > 0
                ? math.asin(absW_clamped)
                : math.PI_DBL - math.asin(absW_clamped);

            // Copy basis vectors
            d.P1[0] = O1.P0;
            d.P1[1] = O1.P1;
            d.P1[2] = O1.P2;
            d.P2[0] = O2.P0;
            d.P2[1] = O2.P1;
            d.P2[2] = O2.P2;
            d.Q1[0] = O1.Q0;
            d.Q1[1] = O1.Q1;
            d.Q1[2] = O1.Q2;
            d.Q2[0] = O2.Q0;
            d.Q2[1] = O2.Q1;
            d.Q2[2] = O2.Q2;

            // Dot products (S vectors = Q * eta)
            d.Pp = Dot(d.P1, d.P2);
            d.Ps = Dot(d.P1, d.Q2) * eta2;
            d.Sp = Dot(d.Q1, d.P2) * eta1;
            d.Ss = Dot(d.Q1, d.Q2) * eta1 * eta2;

            d.P1w = d.P1[0] * w0 + d.P1[1] * w1 + d.P1[2] * w2;
            d.P2w = d.P2[0] * w0 + d.P2[1] * w1 + d.P2[2] * w2;
            d.Q1w = d.Q1[0] * w0 + d.Q1[1] * w1 + d.Q1[2] * w2;
            d.Q2w = d.Q2[0] * w0 + d.Q2[1] * w1 + d.Q2[2] * w2;

            return d;
        }

        private static double Dot(double* a, double* b)
            => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    }
}
