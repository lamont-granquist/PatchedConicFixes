using System.Runtime.CompilerServices;
using Unity.Burst;
using Unity.Mathematics;
using Random = Unity.Mathematics.Random;

// ReSharper disable CompareOfFloatsByEqualityOperator

namespace PatchedConicFixes
{
    [BurstCompile]
    public unsafe struct Baluev
    {
        private const double EPS = 2.2204460492503131e-16;

        /// <summary>
        ///     Core g(u) evaluation shared by all variants.
        ///     The sign parameters encode the ee/he/eh/hh differences:
        ///     signA2: +1 for elliptic second orbit, -1 for hyperbolic second orbit
        ///     signM:  controls sign flips on Sp/Ss terms in M, N for hyperbolic first orbit
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double FuncG(double A, double B, double C, double M, double N,
            double K, double signA2)
        {
            double C2 = C * C;
            double A2 = signA2 * A * A;
            double B2 = B * B;
            double Ac = A2 - C2;
            double Bc = B2 - C2;

            double sNA = signA2; // sign applied to N*A terms mirrors signA2

            return K * K * (Ac * Bc)
                + 2.0 * K * C * (sNA * N * A * Ac + M * B * Bc)
                - (A2 + B2) * (sNA * N * N * Ac + M * M * Bc - 2.0 * sNA * N * M * A * B);
        }

        /// <summary>
        ///     Evaluates g(u) at the j-th sample point for the EE case (both orbits elliptic).
        ///     Returns a real value.
        /// </summary>
        private static double FuncG_EE_Fast(int j, ref SAuxData data, ref TrigTable trig)
        {
            double x = trig.Cos[j];
            double y = trig.Sin[j];

            double A  = data.Ps * y - data.Ss * x;
            double B  = data.Pp * y - data.Sp * x;
            double C  = data.e2 * B - data.alpha1 * data.e1 * y * (1.0 - data.e1 * x);
            double x_ = x - data.e1;
            double M  = data.Sp * y + data.Pp * x_ + data.alpha2 * data.e2;
            double N  = -data.Ss * y - data.Ps * x_;

            return FuncG(A, B, C, M, N, data.K, 1.0);
        }

        /// <summary>
        ///     Evaluates g(u) at the j-th sample point for the EH case
        ///     (first orbit elliptic, second hyperbolic).
        ///     Returns a real value.
        /// </summary>
        private static double FuncG_EH_Fast(int j, ref SAuxData data, ref TrigTable trig)
        {
            double x = trig.Cos[j];
            double y = trig.Sin[j];

            double A  = data.Ps * y - data.Ss * x;
            double B  = data.Pp * y - data.Sp * x;
            double C  = data.e2 * B - data.alpha1 * data.e1 * y * (1.0 - data.e1 * x);
            double x_ = x - data.e1;
            double M  = data.Sp * y + data.Pp * x_ + data.alpha2 * data.e2;
            double N  = -data.Ss * y - data.Ps * x_;

            return FuncG(A, B, C, M, N, data.K, -1.0);
        }

        /// <summary>
        ///     Evaluates g(u) at the j-th sample point for the HE case
        ///     (first orbit hyperbolic, second elliptic).
        ///     Returns a complex value (the hyperbolic substitution u→iu makes intermediates complex).
        /// </summary>
        private static Cmplx FuncG_HE_Fast(int j, ref SAuxData data, ref TrigTable trig)
        {
            double x = trig.Cos[j]; // cos of the imaginary argument
            double y = trig.Sin[j]; // sin of the imaginary argument

            // Hyperbolic first orbit: sin(iu) = i·sinh(u), cos(iu) = cosh(u)
            // So A, B, M, N become complex with real and imaginary parts split
            var    A  = new Cmplx(data.Ps * y, -data.Ss * x);
            var    B  = new Cmplx(data.Pp * y, -data.Sp * x);
            Cmplx  C  = data.e2 * B - data.alpha1 * data.e1 * y * (1.0 - data.e1 * x);
            double x_ = x - data.e1;
            var    M  = new Cmplx(data.Pp * x_ + data.alpha2 * data.e2, data.Sp * y);
            var    N  = new Cmplx(-data.Ps * x_, -data.Ss * y);

            Cmplx C2 = C * C;
            Cmplx A2 = A * A;
            Cmplx B2 = B * B;
            Cmplx Ac = A2 - C2;
            Cmplx Bc = B2 - C2;

            return data.K * data.K * (Ac * Bc)
                + 2.0 * data.K * C * (N * A * Ac + M * B * Bc)
                - (A2 + B2) * (N * N * Ac + M * M * Bc - 2.0 * N * M * A * B);
        }

        /// <summary>
        ///     Evaluates g(u) at the j-th sample point for the HH case
        ///     (both orbits hyperbolic).
        ///     Returns a complex value.
        /// </summary>
        private static Cmplx FuncG_HH_Fast(int j, ref SAuxData data, ref TrigTable trig)
        {
            double x = trig.Cos[j];
            double y = trig.Sin[j];

            var    A  = new Cmplx(data.Ps * y, -data.Ss * x);
            var    B  = new Cmplx(data.Pp * y, -data.Sp * x);
            Cmplx  C  = data.e2 * B - data.alpha1 * data.e1 * y * (1.0 - data.e1 * x);
            double x_ = x - data.e1;
            var    M  = new Cmplx(data.Pp * x_ + data.alpha2 * data.e2, data.Sp * y);
            var    N  = new Cmplx(-data.Ps * x_, -data.Ss * y);

            Cmplx C2 = C * C;
            Cmplx A2 = -(A * A); // negated for hyperbolic second orbit
            Cmplx B2 = B * B;
            Cmplx Ac = A2 - C2;
            Cmplx Bc = B2 - C2;

            return data.K * data.K * (Ac * Bc)
                + 2.0 * data.K * C * (-N * A * Ac + M * B * Bc)
                - (A2 + B2) * (-N * N * Ac + M * M * Bc + 2.0 * N * M * A * B);
        }

        private const int DEG = 16;

        /// <summary>
        ///     Builds the degree-16 algebraic polynomial from the trigonometric polynomial g(u)
        ///     via DFT sampling at 21 points.
        ///     For elliptic first orbit (EE or EH): g(u_j) is real, coefficients are complex
        ///     with conjugation symmetry c[k] = conj(c[DEG-k]).
        ///     For hyperbolic first orbit (HE or HH): g(u_j) is complex, coefficients are real.
        ///     Returns the estimated coefficient uncertainty (cerr).
        /// </summary>
        private static double CreatePolynomial(ref SAuxData data, ref TrigTable trig,
            Cmplx* c, OrbitPairType type)
        {
            if (type == OrbitPairType.EE || type == OrbitPairType.EH)
                return CreatePolynomialElliptic(ref data, ref trig, c, type);
            return CreatePolynomialHyperbolic(ref data, ref trig, c, type);
        }

        /// <summary>
        ///     Elliptic first orbit (e1 ≤ 1). g(u) is real-valued.
        ///     Polynomial coefficients are complex with c[k] = conj(c[DEG-k]).
        /// </summary>
        private static double CreatePolynomialElliptic(ref SAuxData data, ref TrigTable trig,
            Cmplx* c, OrbitPairType type)
        {
            double ee    = data.e1 * data.e2;
            double scale = data.alpha1 * data.e1 * data.e1 / 16.0;

            // Direct evaluation of the highest coefficient c[0] = hc, c[DEG] = conj(hc)
            Cmplx hc;
            if (type == OrbitPairType.EE)
            {
                hc = scale * scale
                    * new Cmplx(data.Pp - data.Ss - ee, data.Sp + data.Ps)
                    * new Cmplx(data.Pp - data.Ss + ee, data.Sp + data.Ps)
                    * new Cmplx(data.Pp + data.Ss - ee, data.Sp - data.Ps)
                    * new Cmplx(data.Pp + data.Ss + ee, data.Sp - data.Ps);
            }
            else // EH
            {
                hc = scale * scale
                    * new Cmplx(data.Pp - data.Ps - ee, data.Sp - data.Ss)
                    * new Cmplx(data.Pp - data.Ps + ee, data.Sp - data.Ss)
                    * new Cmplx(data.Pp + data.Ps - ee, data.Sp + data.Ss)
                    * new Cmplx(data.Pp + data.Ps + ee, data.Sp + data.Ss);
            }

            // Sample g(u) at 21 points
            double* vals = stackalloc double[TrigTable.DIM];
            for (int j = 0; j < TrigTable.DIM; j++)
            {
                vals[j] = type == OrbitPairType.EE
                    ? FuncG_EE_Fast(j, ref data, ref trig)
                    : FuncG_EH_Fast(j, ref data, ref trig);
            }

            // Preprocess: split into even (sum) and odd (difference) parts
            for (int j = 1; j <= TrigTable.DIM / 2; j++)
            {
                double sum = vals[j] + vals[TrigTable.DIM - j];
                double dif = vals[j] - vals[TrigTable.DIM - j];
                vals[j]                 = sum;
                vals[TrigTable.DIM - j] = dif;
            }

            // DFT for the middle coefficient c[DEG/2] (k=0 term)
            double a = vals[0];
            for (int j = 1; j <= TrigTable.DIM / 2; j++)
                a += vals[j];
            c[DEG / 2] = a / TrigTable.DIM;

            // DFT for coefficients c[DEG/2-i] and c[DEG/2+i], i = 1..7
            int i;
            for (i = 1; i < DEG / 2; i++)
            {
                a = vals[0];
                double b = 0.0;
                for (int j = 1; j <= TrigTable.DIM / 2; j++)
                {
                    int idx = j * i % TrigTable.DIM;
                    a += vals[j] * trig.Cos[idx];
                    b += vals[TrigTable.DIM - j] * trig.Sin[idx];
                }

                a              /= TrigTable.DIM;
                b              /= TrigTable.DIM;
                c[DEG / 2 - i] =  new Cmplx(a, b);
                c[DEG / 2 + i] =  new Cmplx(a, -b);
            }

            // Set highest/lowest coefficients from direct formula
            c[0]   = hc;
            c[DEG] = hc.Conjugate();

            // Error estimation: compare DFT estimate of c[8] against the analytic hc
            double cerr = 0.0;

            // i == DEG/2 here (continuing from the loop)
            a = vals[0];
            double b2 = 0.0;
            for (int j = 1; j <= TrigTable.DIM / 2; j++)
            {
                int idx = j * i % TrigTable.DIM;
                a  += vals[j] * trig.Cos[idx];
                b2 += vals[TrigTable.DIM - j] * trig.Sin[idx];
            }

            a    /= TrigTable.DIM;
            b2   /= TrigTable.DIM;
            cerr += (new Cmplx(a, b2) - hc).Norm;

            // Remaining aliased coefficients (should be zero; their magnitude estimates error)
            for (i++; i <= TrigTable.DIM / 2; i++)
            {
                a = vals[0];
                double b = 0.0;
                for (int j = 1; j <= TrigTable.DIM / 2; j++)
                {
                    int idx = j * i % TrigTable.DIM;
                    a += vals[j] * trig.Cos[idx];
                    b += vals[TrigTable.DIM - j] * trig.Sin[idx];
                }

                a    /= TrigTable.DIM;
                b    /= TrigTable.DIM;
                cerr += a * a + b * b;
            }

            cerr /= (TrigTable.DIM - DEG + 1) / 2.0;

            return math.max(math.sqrt(cerr), EPS * hc.Abs);
        }

        public static double Sqr(double x) => x * x;

        /// <summary>
        ///     Hyperbolic first orbit (e1 > 1). g(u) is complex-valued.
        ///     Polynomial coefficients are real.
        /// </summary>
        private static double CreatePolynomialHyperbolic(ref SAuxData data, ref TrigTable trig,
            Cmplx* c, OrbitPairType type)
        {
            double d   = data.e2 * data.e2;
            double d1p = Sqr((data.Pp + data.Sp) / data.e1);
            double d2p = Sqr((data.Ps + data.Ss) / data.e1);
            double dpp = d1p + d2p;
            double dpm = d1p - d2p;

            double prefactor = Sqr(data.alpha1 * Sqr(Sqr(data.e1)) / 16.0);

            double h = prefactor * (type == OrbitPairType.HE
                ? dpp * dpp + d * d - 2.0 * d * dpm
                : dpm * dpm + d * d - 2.0 * d * dpp);

            double d1m = Sqr((data.Pp - data.Sp) / data.e1);
            double d2m = Sqr((data.Ps - data.Ss) / data.e1);
            double dmp = d1m + d2m;
            double dmm = d1m - d2m;

            double hc = prefactor * (type == OrbitPairType.HE
                ? dmp * dmp + d * d - 2.0 * d * dmm
                : dmm * dmm + d * d - 2.0 * d * dmp);

            // Sample g(u) at DIM/2+1 points (symmetry: second half is conjugate of first in reverse)
            Cmplx* vals = stackalloc Cmplx[TrigTable.DIM / 2 + 1];
            double vals0 = (type == OrbitPairType.HE
                ? FuncG_HE_Fast(0, ref data, ref trig)
                : FuncG_HH_Fast(0, ref data, ref trig)).Re; // g(0) is always real

            for (int j = 1; j <= TrigTable.DIM / 2; j++)
            {
                vals[j] = type == OrbitPairType.HE
                    ? FuncG_HE_Fast(j, ref data, ref trig)
                    : FuncG_HH_Fast(j, ref data, ref trig);
            }

            // DFT for the middle coefficient (k=0)
            double a = 0.0;
            for (int j = 1; j <= TrigTable.DIM / 2; j++)
                a += vals[j].Re;
            c[DEG / 2] = (2.0 * a + vals0) / TrigTable.DIM;

            // DFT for coefficients c[DEG/2±i], i = 1..7
            int i;
            for (i = 1; i < DEG / 2; i++)
            {
                a = 0.0;
                double b = 0.0;
                for (int j = 1; j <= TrigTable.DIM / 2; j++)
                {
                    int idx = j * i % TrigTable.DIM;
                    a += vals[j].Re * trig.Cos[idx];
                    b += vals[j].Im * trig.Sin[idx];
                }

                c[DEG / 2 - i] = (2.0 * (a - b) + vals0) / TrigTable.DIM;
                c[DEG / 2 + i] = (2.0 * (a + b) + vals0) / TrigTable.DIM;
            }

            // Set highest/lowest from direct formula
            c[0]   = hc;
            c[DEG] = h;

            // Error estimation
            a = 0.0;
            double b2 = 0.0;
            for (int j = 1; j <= TrigTable.DIM / 2; j++)
            {
                int idx = j * i % TrigTable.DIM;
                a  += vals[j].Re * trig.Cos[idx];
                b2 += vals[j].Im * trig.Sin[idx];
            }

            double cerr = Sqr((2.0 * (a - b2) + vals0) / TrigTable.DIM - hc)
                + Sqr((2.0 * (a + b2) + vals0) / TrigTable.DIM - h);

            for (i++; i <= TrigTable.DIM / 2; i++)
            {
                a = 0.0;
                double b = 0.0;
                for (int j = 1; j <= TrigTable.DIM / 2; j++)
                {
                    int idx = j * i % TrigTable.DIM;
                    a += vals[j].Re * trig.Cos[idx];
                    b += vals[j].Im * trig.Sin[idx];
                }

                cerr += Sqr((2.0 * a + vals0) / TrigTable.DIM) + Sqr(2.0 * b / TrigTable.DIM);
            }

            cerr /= TrigTable.DIM - DEG + 1;

            return math.max(math.sqrt(cerr), EPS * math.max(math.abs(hc), math.abs(h)));
        }

        private enum OrbitPairType
        {
            EE, // both elliptic
            EH, // first elliptic, second hyperbolic
            HE, // first hyperbolic, second elliptic
            HH  // both hyperbolic
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double AngleWrap(double x)
        {
            const double TwoPi = 2.0 * math.PI_DBL;
            x = (x + math.PI_DBL) % TwoPi;
            if (x < 0) x += TwoPi;
            return x - math.PI_DBL;
        }

        /// <summary>
        ///     Hyperbolic atan2: returns the real part of atanh(y/x) = (log|x+y| - log|x-y|) / 2.
        ///     Used for the hyperbolic eccentric anomaly where the circular atan2 is replaced
        ///     by its hyperbolic counterpart.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Atan2h(double y, double x) => (math.log(math.abs(x + y)) - math.log(math.abs(x - y))) / 2.0;

        /// <summary>
        ///     Returns atan2(y, x) for elliptic orbits (trig=true)
        ///     or atan2h(y, x) for hyperbolic orbits (trig=false).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Atan2Smart(double y, double x, bool trig) => trig ? math.atan2(y, x) : Atan2h(y, x);

        /// <summary>
        ///     Computes the position vector on an orbit at eccentric anomaly u.
        ///     r = a * [(cos u - e)·P + sin u · η · Q]          (elliptic)
        ///     r = a * [(cosh u - e)·P + (-sinh u) · η · Q]     (hyperbolic)
        ///     The 'br' flag selects the branch for hyperbolic orbits.
        ///     For elliptic orbits br is always true.
        ///     When br=false, the position is negated (opposite branch).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void RadiusVector(bool br, double* P, double* Q,
            double a, double e, double u, double* r)
        {
            double p, q;
            if (e <= 1.0)
            {
                p = math.cos(u);
                q = math.sin(u) * math.sqrt(1.0 - e * e);
            }
            else
            {
                double tmp = math.exp(math.abs(u));
                p = (tmp + 1.0 / tmp) / 2.0;                                         // cosh(|u|)
                q = (1.0 / tmp - tmp) / 2.0 * math.sign(u) * math.sqrt(e * e - 1.0); // -sinh(u)·η
            }

            if (!br)
            {
                p = -p;
                q = -q;
            }

            p    -= e;
            r[0] =  a * (P[0] * p + Q[0] * q);
            r[1] =  a * (P[1] * p + Q[1] * q);
            r[2] =  a * (P[2] * p + Q[2] * q);
        }

        /// <summary>
        ///     Computes the position vector r and the derivative dr/du (divided by a).
        ///     rd = dr/du / a = -sin u · P + cos u · η · Q      (elliptic)
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void RadiusVectorValDer(bool br, double* P, double* Q,
            double a, double e, double u,
            double* r, double* rd)
        {
            double x, y, eta, h;
            if (e <= 1.0)
            {
                x   = math.cos(u);
                y   = math.sin(u);
                eta = math.sqrt(1.0 - e * e);
                h   = x * eta;
            }
            else
            {
                double tmp = math.exp(math.abs(u));
                x   = (tmp + 1.0 / tmp) / 2.0;
                y   = (1.0 / tmp - tmp) / 2.0 * math.sign(u);
                eta = math.sqrt(e * e - 1.0);
                h   = -x * eta;
            }

            if (!br)
            {
                x = -x;
                y = -y;
                h = -h;
            }

            double p = x - e;
            double q = y * eta;

            r[0] = a * (P[0] * p + Q[0] * q);
            r[1] = a * (P[1] * p + Q[1] * q);
            r[2] = a * (P[2] * p + Q[2] * q);

            rd[0] = -P[0] * y + Q[0] * h;
            rd[1] = -P[1] * y + Q[1] * h;
            rd[2] = -P[2] * y + Q[2] * h;
        }

        /// <summary>
        ///     Computes the position vector r, derivative dr/du / a, and second derivative d²r/du² / a.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static void RadiusVectorValDer2(bool br, double* P, double* Q,
            double a, double e, double u,
            double* r, double* rd, double* rdd)
        {
            double x, y, eta, h;
            if (e <= 1.0)
            {
                x   = math.cos(u);
                y   = math.sin(u);
                eta = math.sqrt(1.0 - e * e);
                h   = x * eta;
            }
            else
            {
                double tmp = math.exp(math.abs(u));
                x   = (tmp + 1.0 / tmp) / 2.0;
                y   = (1.0 / tmp - tmp) / 2.0 * math.sign(u);
                eta = math.sqrt(e * e - 1.0);
                h   = -x * eta;
            }

            if (!br)
            {
                x = -x;
                y = -y;
                h = -h;
            }

            double p = x - e;
            double q = y * eta;

            r[0] = a * (P[0] * p + Q[0] * q);
            r[1] = a * (P[1] * p + Q[1] * q);
            r[2] = a * (P[2] * p + Q[2] * q);

            rd[0] = -P[0] * y + Q[0] * h;
            rd[1] = -P[1] * y + Q[1] * h;
            rd[2] = -P[2] * y + Q[2] * h;

            rdd[0] = P[0] * x + Q[0] * q;
            rdd[1] = P[1] * x + Q[1] * q;
            rdd[2] = P[2] * x + Q[2] * q;
            if (e <= 1.0)
            {
                rdd[0] = -rdd[0];
                rdd[1] = -rdd[1];
                rdd[2] = -rdd[2];
            }
        }

        /// <summary>
        ///     Computes the Euclidean distance between two points on the orbits
        ///     at eccentric anomalies u1 and u2.
        /// </summary>
        public static double DistanceBetween(SAuxData* data, bool br1, bool br2,
            double u1, double u2)
        {
            double* r1 = stackalloc double[3];
            double* r2 = stackalloc double[3];

            RadiusVector(br1, data->P1, data->Q1, data->a1, data->e1, u1, r1);
            RadiusVector(br2, data->P2, data->Q2, data->a2, data->e2, u2, r2);

            double dx = r2[0] - r1[0];
            double dy = r2[1] - r1[1];
            double dz = r2[2] - r1[2];

            return math.sqrt(dx * dx + dy * dy + dz * dz);
        }

        /// <summary>
        ///     Analytically solves for the eccentric anomaly u2 on the second orbit,
        ///     given u1 on the first orbit, at a critical point of the distance function.
        ///     Produces two candidate values u2 and u2Alt. At a true critical point,
        ///     one of them (most likely u2) corresponds to the actual critical point.
        ///     The gradient condition M·sin(u2) + N·cos(u2) = K·sin(u2)·cos(u2) is used
        ///     to select which candidate is the better match.
        ///     Returns true if the discriminant was negative (degenerate case — both outputs
        ///     are set to the same fallback value). Returns false on success.
        ///     If perror is non-null, it is updated with the propagated error estimate.
        /// </summary>
        public static bool EliminatedAnomaly(SAuxData* data, double u1,
            out double u2, out double u2Alt,
            double* perror = null)
        {
            double x, y, x_, M, N;

            if (data->e1 <= 1.0)
            {
                x  = math.cos(u1);
                y  = math.sin(u1);
                x_ = x - data->e1;
                M  = data->Sp * y + data->Pp * x_ + data->alpha2 * data->e2;
                N  = -data->Ss * y - data->Ps * x_;
            }
            else
            {
                double tmp1 = math.exp(math.abs(u1));
                x  = (tmp1 + 1.0 / tmp1) / 2.0;
                y  = (tmp1 - 1.0 / tmp1) / 2.0 * math.sign(u1);
                x_ = x - data->e1;
                M  = -data->Sp * y + data->Pp * x_ + data->alpha2 * data->e2;
                N  = data->Ss * y - data->Ps * x_;
            }

            double A = data->Ps * y - data->Ss * x;
            double B = data->Pp * y - data->Sp * x;
            double C = data->e2 * B - data->alpha1 * data->e1 * y * (1.0 - data->e1 * x);

            double V, D;
            if (data->e2 <= 1.0)
            {
                V = A * A + B * B;
                D = V - C * C;
            }
            else
            {
                V = B * B - A * A;
                D = C * C - V;
            }

            bool negativeDiscriminant = D < 0;
            if (negativeDiscriminant)
            {
                D = 0;
                if (data->e2 <= 1.0)
                    u2 = C > 0 ? math.atan2(A, B) : math.atan2(-A, -B);
                else
                    u2 = Atan2h(A, B);
                u2Alt = u2;
            }
            else
            {
                D = math.sqrt(D);
                double AC = A * C;
                double BD = B * D;
                double y1 = (AC - BD) / V;
                double y2 = (AC + BD) / V;
                if (!(data->e2 <= 1.0))
                {
                    double t = y1;
                    y1 = y2;
                    y2 = t;
                }

                double BC = B * C;
                double AD = A * D;
                double x1 = (BC + AD) / V;
                double x2 = (BC - AD) / V;

                if (data->e2 <= 1.0)
                {
                    u2    = math.atan2(y1, x1);
                    u2Alt = math.atan2(y2, x2);
                }
                else
                {
                    if (!(x1 > 0) && x2 > 0)
                    {
                        double t;
                        t  = x1;
                        x1 = x2;
                        x2 = t;
                        t  = y1;
                        y1 = y2;
                        y2 = t;
                    }

                    u2    = Atan2h(y1, x1);
                    u2Alt = Atan2h(y2, x2);
                }

                // Select which candidate better satisfies the gradient condition
                double err  = M * y1 + N * x1 - data->K * y1 * x1;
                double err2 = M * y2 + N * x2 - data->K * y2 * x2;
                if (!(math.abs(err) < math.abs(err2)) && (data->e2 <= 1.0 || x2 > 0))
                {
                    double t = u2;
                    u2    = u2Alt;
                    u2Alt = t;
                }
            }

            if (perror != null)
            {
                double sigma = *perror * math.sqrt(1.0 + Sqr(data->e2) + Sqr(data->alpha1 * data->e1));
                *perror = sigma / math.sqrt(D * D + math.abs(C) * sigma / 2.0);
            }

            return negativeDiscriminant;
        }

        /// <summary>
        ///     Computes the candidate minimum distance for a given eccentric anomaly u1
        ///     on the first orbit. The second anomaly u2 is determined analytically via
        ///     EliminatedAnomaly, and the smaller of the two candidate distances is returned.
        ///     Returns a negative value if the discriminant was negative (degenerate case).
        ///     If perror is non-null, it is updated with the propagated error estimate.
        /// </summary>
        public static double MinDistanceFor(SAuxData* data, double u1, out double u2,
            double* perror = null)
        {
            double u2a, u2b;
            if (EliminatedAnomaly(data, u1, out u2a, out u2b, perror))
            {
                u2 = u2a;
                return -DistanceBetween(data, true, true, u1, u2a);
            }

            double da = DistanceBetween(data, true, true, u1, u2a);
            double db = DistanceBetween(data, true, true, u1, u2b);

            if (da < db)
            {
                u2 = u2a;
                return da;
            }

            u2 = u2b;
            return db;
        }

        // ============================================================
        // Layer 4: 2D squared-distance surface
        // ============================================================

        /// <summary>
        ///     Computes the normalized squared distance ρ(u1,u2) and its gradient.
        ///     ρ = |r2 - r1|² / (2·a1·a2)
        ///     g[0] = ∂ρ/∂u1
        ///     g[1] = ∂ρ/∂u2
        ///     The normalization by 2·a1·a2 keeps the intermediate quantities
        ///     from scaling with the orbital radii.
        /// </summary>
        public static double SqDistVG(SAuxData* data, bool br1, bool br2,
            double u1, double u2, double* g)
        {
            double* r1  = stackalloc double[3];
            double* rd1 = stackalloc double[3];
            double* r2  = stackalloc double[3];
            double* rd2 = stackalloc double[3];

            RadiusVectorValDer(br1, data->P1, data->Q1, data->a1, data->e1, u1, r1, rd1);
            RadiusVectorValDer(br2, data->P2, data->Q2, data->a2, data->e2, u2, r2, rd2);

            double dr0 = r2[0] - r1[0];
            double dr1 = r2[1] - r1[1];
            double dr2 = r2[2] - r1[2];

            g[0] = -(dr0 * rd1[0] + dr1 * rd1[1] + dr2 * rd1[2]) / data->a2;
            g[1] = (dr0 * rd2[0] + dr1 * rd2[1] + dr2 * rd2[2]) / data->a1;

            return (dr0 * dr0 + dr1 * dr1 + dr2 * dr2) / (2.0 * data->a1 * data->a2);
        }

        /// <summary>
        ///     Computes the normalized squared distance ρ(u1,u2), its gradient, and its Hessian.
        ///     H[0] = ∂²ρ/∂u1²
        ///     H[1] = ∂²ρ/∂u2²
        ///     H[2] = ∂²ρ/∂u1∂u2
        ///     The Hessian is stored as three independent elements of the symmetric 2×2 matrix.
        /// </summary>
        public static double SqDistVGH(SAuxData* data, bool br1, bool br2,
            double u1, double u2, double* g, double* H)
        {
            double* r1   = stackalloc double[3];
            double* rd1  = stackalloc double[3];
            double* rdd1 = stackalloc double[3];
            double* r2   = stackalloc double[3];
            double* rd2  = stackalloc double[3];
            double* rdd2 = stackalloc double[3];

            RadiusVectorValDer2(br1, data->P1, data->Q1, data->a1, data->e1, u1, r1, rd1, rdd1);
            RadiusVectorValDer2(br2, data->P2, data->Q2, data->a2, data->e2, u2, r2, rd2, rdd2);

            double dr0 = r2[0] - r1[0];
            double dr1 = r2[1] - r1[1];
            double dr2 = r2[2] - r1[2];

            double dotDrRd1  = dr0 * rd1[0] + dr1 * rd1[1] + dr2 * rd1[2];
            double dotDrRd2  = dr0 * rd2[0] + dr1 * rd2[1] + dr2 * rd2[2];
            double dotDrRdd1 = dr0 * rdd1[0] + dr1 * rdd1[1] + dr2 * rdd1[2];
            double dotDrRdd2 = dr0 * rdd2[0] + dr1 * rdd2[1] + dr2 * rdd2[2];
            double normRd1   = rd1[0] * rd1[0] + rd1[1] * rd1[1] + rd1[2] * rd1[2];
            double normRd2   = rd2[0] * rd2[0] + rd2[1] * rd2[1] + rd2[2] * rd2[2];
            double dotRd1Rd2 = rd1[0] * rd2[0] + rd1[1] * rd2[1] + rd1[2] * rd2[2];

            g[0] = -dotDrRd1 / data->a2;
            g[1] = dotDrRd2 / data->a1;

            H[0] = (-dotDrRdd1 + normRd1 * data->a1) / data->a2;
            H[1] = (dotDrRdd2 + normRd2 * data->a2) / data->a1;
            H[2] = -dotRd1Rd2;

            return (dr0 * dr0 + dr1 * dr1 + dr2 * dr2) / (2.0 * data->a1 * data->a2);
        }

        /// <summary>
        ///     Performs a single 2D Newton step on the squared distance surface ρ(u1,u2).
        ///     Computes ρ, gradient, Hessian, then solves H·δu = -g and updates u1, u2.
        ///     Outputs:
        ///     rho  - value of ρ before the step
        ///     g    - gradient before the step
        ///     H    - Hessian before the step
        ///     detH - determinant of the Hessian
        ///     Returns |δu|² = δu1² + δu2² (the squared step size), or 0 if detH = 0.
        /// </summary>
        public static double SqDist2DIter(SAuxData* data, bool br1, bool br2,
            ref double u1, ref double u2,
            out double rho, double* g, double* H, out double detH)
        {
            rho  = SqDistVGH(data, br1, br2, u1, u2, g, H);
            detH = H[0] * H[1] - H[2] * H[2];
            if (detH == 0.0) return 0.0;

            double du1 = -(H[1] * g[0] - H[2] * g[1]) / detH;
            double du2 = -(H[0] * g[1] - H[2] * g[0]) / detH;

            u1 += du1;
            u2 += du2;

            return du1 * du1 + du2 * du2;
        }


        /// <summary>
        ///     Estimates the numerical uncertainty of a root of the polynomial.
        ///     Uses a quadratic approximation based on the polynomial value, first, and second
        ///     derivatives at the root, plus the propagated coefficient error cerr.
        ///     The quadratic approximation solves f(z+δ) ≈ f + f'δ + f''δ²/2 = 0 for |δ|,
        ///     using the discriminant D = sqrt(f'² - 2f·f'') to pick the more numerically
        ///     stable root of the quadratic.
        ///     Returns the estimated relative uncertainty of the root.
        /// </summary>
        public static double RootError(int n, Cmplx* c, double cerr, Cmplx z)
        {
            double m        = z.Norm; // |z|²
            bool   forward  = m <= 1.0;
            if (!forward) m = 1.0 / m;

            // Expected error in polynomial value from coefficient errors
            // ferr² = cerr² · Σ_{i=0}^{n} |z|^{2i} = cerr² · (1 + m + m² + ... + m^n)
            double ferr2 = 1.0;
            for (int i = 1; i <= n; i++)
            {
                ferr2 *= m;
                ferr2 += 1.0;
            }

            ferr2 *= cerr * cerr;

            // Evaluate polynomial, first, and second derivatives
            PolySolver.PolynomialValDer2(n, c, forward ? z : z.Inverse(), out Cmplx f, out Cmplx fd, out Cmplx fd2, forward);

            // Quadratic approximation: pick the discriminant branch with larger denominator
            var D = Cmplx.Sqrt(fd * fd - 2.0 * f * fd2);

            Cmplx  fdPlusD  = fd + D;
            Cmplx  fdMinusD = fd - D;
            double d2       = 4.0 * f.Norm / math.max(fdPlusD.Norm, fdMinusD.Norm);

            return math.sqrt((d2 + ferr2 / D.Norm) / m);
        }

        /// <summary>
        ///     Estimates an upper bound on the absolute value of all roots of the polynomial.
        ///     Uses the Cauchy bound: max over i of (|c[n-i]/c[n]|)^{1/i} for forward,
        ///     or the equivalent for backward.
        ///     Forward mode gives an upper bound on |z| for roots.
        ///     Backward mode gives an upper bound on |1/z|, i.e. a lower bound on |z|.
        ///     Returns the bound (as a modulus, not squared).
        /// </summary>
        public static double MaxRootBound(int n, Cmplx* c, bool forward)
        {
            double nch = forward ? c[n].Norm : c[0].Norm;
            if (n < 1) return double.PositiveInfinity;

            double res = (forward ? c[n - 1].Norm : c[1].Norm) / nch;

            for (int i = 2; i <= n; i++)
            {
                double ci          = forward ? c[n - i].Norm : c[i].Norm;
                double tmp         = math.pow(ci / nch, 1.0 / i);
                if (res < tmp) res = tmp;
            }

            return math.sqrt(res);
        }

        /// <summary>
        ///     Refines the critical point (u1, u2) of the squared distance surface ρ(u1, u2)
        ///     using 2D Newton iteration, and estimates residual uncertainties.
        ///     Iterates until the step size stops improving or the desired precision eps is reached.
        ///     Outputs:
        ///     u1, u2      - refined eccentric anomalies
        ///     u1Err, u2Err - estimated numeric uncertainties in the anomalies
        ///     rho         - final value of ρ
        ///     rhoErr      - estimated numeric uncertainty in ρ
        ///     iterCnt     - number of iterations performed
        ///     g           - final gradient (2 elements)
        ///     H           - final Hessian (3 elements: H00, H11, H01)
        ///     Returns Hessian definiteness indicator:
        ///     +2 if positive definite (true minimum)
        ///     +1 if positive semi-definite
        ///     -2 if negative definite (maximum)
        ///     -1 if negative semi-definite
        ///     0 if indefinite (saddle point)
        /// </summary>
        public static short NewtonSqDist(SAuxData* data,
            ref double u1, ref double u2,
            out double u1Err, out double u2Err,
            out double rho, out double rhoErr,
            double eps, out ulong iterCnt, int maxcount,
            double* g, double* H)
        {
            double detH;
            double derr2;
            double err2 = double.PositiveInfinity;

            iterCnt = 0;
            do
            {
                double err2_ = SqDist2DIter(data, true, true, ref u1, ref u2,
                    out rho, g, H, out detH);
                if (data->e1 <= 1.0 && math.abs(u1) > math.PI_DBL) u1 = AngleWrap(u1);
                if (data->e2 <= 1.0 && math.abs(u2) > math.PI_DBL) u2 = AngleWrap(u2);
                if (detH == 0.0) break;
                derr2 = err2 - err2_;
                err2  = err2_;
                iterCnt++;
            } while (derr2 > EPS * err2 && err2 > eps * eps && iterCnt <= (ulong)maxcount);

            // Final evaluation for residual gradient
            rho = SqDistVG(data, true, true, u1, u2, g);

            // Error estimates from residual gradient and Hessian
            u1Err = math.abs((H[1] * g[0] - H[2] * g[1]) / detH);
            u2Err = math.abs((H[0] * g[1] - H[2] * g[0]) / detH);
            rhoErr = math.abs((H[1] * Sqr(g[0]) + H[0] * Sqr(g[1])
                - 2.0 * H[2] * g[0] * g[1]) / (2.0 * detH));

            // Degenerate Hessian
            if (detH == 0.0)
            {
                u1Err  = 0;
                u2Err  = 0;
                rhoErr = 0;
                if (H[0] > 0 || H[1] > 0) return +1;
                if (H[0] < 0 || H[1] < 0) return -1;
            }

            // Classify definiteness
            if (detH > 0)
            {
                if (H[0] > 0) return +2;
                if (H[0] < 0) return -2;
            }

            return 0;
        }

        // this is a bit tightly coupled to CheckEncounter, but so is FindAllMinima()
        public struct MoidInfo
        {
            public double dst;
            public double u1;
            public double u2;
            public double ta1;
            public double ta2;
            public double tt;
            public double ut;
        }

        public static int FindAllMinima(in COrbitData O1, in COrbitData O2, MoidInfo* info)
        {
            const double MAX_ROOT_ERR = 1.49012e-08;
            const double MIN_ROOT_ERR = 4.44089e-16;
            const double NU           = 1.0;

            int rootCount = 0;

            bool       swapped = O1.e < O2.e;
            COrbitData A       = swapped ? O2 : O1;
            COrbitData B       = swapped ? O1 : O2;
            var        data    = SAuxData.Create(in A, in B);

            // can't solve two perfectly circular orbits
            if (!(data.e1 != 1.0 && data.e2 != 1.0))
                return 0;

            var trig = TrigTable.Create();

            Cmplx* c  = stackalloc Cmplx[DEG + 1];
            Cmplx* c_ = stackalloc Cmplx[DEG + 1];

            // Build the polynomial and estimate coefficient error
            OrbitPairType type;
            if (data.e1 <= 1.0)
                type = data.e2 <= 1.0 ? OrbitPairType.EE : OrbitPairType.EH;
            else
                type = data.e2 <= 1.0 ? OrbitPairType.HE : OrbitPairType.HH;

            double cerr = CreatePolynomial(ref data, ref trig, c, type);

            // Rescale coefficients for hyperbolic case (geometric mean normalization)
            double znrm = RescaleCoefficients(data, c, ref cerr);

            // Copy coefficients for root finding (c_ gets destroyed by deflation, c kept for error estimation)
            for (int i = 0; i <= DEG; i++)
                c_[i] = c[i];

            Cmplx* roots = stackalloc Cmplx[DEG];

            var rng = new Random(42);

            GenerateRootGuess(roots, c, data, znrm, rng);

            // --- Find all 16 roots ---
            PolySolver.PolynomialRoots(DEG, c_, roots, MIN_ROOT_ERR, MAX_ROOT_ERR, data.e1 <= 1.0, rng);

            double* g = stackalloc double[2];
            double* H = stackalloc double[3];

            for (int i = 0; i < DEG; i++)
            {
                // we should be absolutely done finding minima after 4, no solutions with 5+ are known.
                if (rootCount >= 4) break;

                double rooterr = RootError(DEG, c, cerr, roots[i]) * NU;

                double dlt;
                if (data.e1 <= 1.0)
                    dlt = math.abs(math.log(roots[i].Norm)) / (2.0 * rooterr);
                else
                    dlt = (roots[i] != 0.0 ? math.abs(roots[i].Phase) : 0.0) / rooterr;

                if (dlt < 3)
                {
                    double u1;
                    if (data.e1 <= 1.0)
                        u1 = roots[i].Phase;
                    else
                        u1 = roots[i].Re > 0 ? -math.log(roots[i].Re * znrm) : double.PositiveInfinity;

                    double d = math.abs(MinDistanceFor(&data, u1, out double u2));
                    if (d < 0) continue; // degenerate root

                    short Hsign = NewtonSqDist(&data, ref u1, ref u2, out _, out _, out double dst, out _, MIN_ROOT_ERR, out _, 30, g, H);
                    Hsign *= (short)(data.a1 * data.a2 >= 0 ? 1 : -1);

                    // 2 == positive definite (minima, not saddle or maxima)
                    if (Hsign != +2) continue;

                    dst = math.sqrt(dst * 2.0 * data.a1 * data.a2);

                    info[rootCount].dst = dst;
                    info[rootCount].u1  = swapped ? u2 : u1;
                    info[rootCount].u2  = swapped ? u1 : u2;

                    rootCount++;
                }
            }

            return rootCount;
        }

        public static SMOIDResult MOID_fast(in COrbitData O1, in COrbitData O2,
            double maxrooterr = 1.49012e-08, double minrooterr = 4.44089e-16,
            double nu = 1.0)
        {
            var result = SMOIDResult.Default();
            var data   = SAuxData.Create(in O1, in O2);

            // can't solve two perfectly circular orbits
            if (!(data.e1 != 1.0 && data.e2 != 1.0))
            {
                result.good = false;
                return result;
            }

            var trig = TrigTable.Create();

            Cmplx* c  = stackalloc Cmplx[DEG + 1];
            Cmplx* c_ = stackalloc Cmplx[DEG + 1];

            // Build the polynomial and estimate coefficient error
            OrbitPairType type;
            if (data.e1 <= 1.0)
                type = data.e2 <= 1.0 ? OrbitPairType.EE : OrbitPairType.EH;
            else
                type = data.e2 <= 1.0 ? OrbitPairType.HE : OrbitPairType.HH;

            double cerr = CreatePolynomial(ref data, ref trig, c, type);

            // Rescale coefficients for hyperbolic case (geometric mean normalization)
            double znrm = RescaleCoefficients(data, c, ref cerr);

            // Copy coefficients for root finding (c_ gets destroyed by deflation, c kept for error estimation)
            for (int i = 0; i <= DEG; i++)
                c_[i] = c[i];

            Cmplx* roots = stackalloc Cmplx[DEG];

            var rng = new Random(42);

            GenerateRootGuess(roots, c, data, znrm, rng);

            // --- Find all 16 roots ---
            result.iterCount = PolySolver.PolynomialRoots(DEG, c_, roots, minrooterr, maxrooterr, data.e1 <= 1.0, rng);

            // --- Filter roots, recover u2, find minimum distance ---
            for (int i = 0; i < DEG; i++)
            {
                double rooterr = RootError(DEG, c, cerr, roots[i]) * nu;

                double dlt;
                if (data.e1 <= 1.0)
                    dlt = math.abs(math.log(roots[i].Norm)) / (2.0 * rooterr);
                else
                    dlt = (roots[i] != 0.0 ? math.abs(roots[i].Phase) : 0.0) / rooterr;

                if (dlt < 10)
                {
                    result.good &= rooterr < maxrooterr;
                    if (dlt < 3)
                    {
                        result.rootCount++;

                        double u1;
                        if (data.e1 <= 1.0)
                            u1 = roots[i].Phase;
                        else
                            u1 = roots[i].Re > 0 ? -math.log(roots[i].Re * znrm) : double.PositiveInfinity;

                        double dst = math.abs(MinDistanceFor(&data, u1, out double u2));

                        if (dst >= 0 && (result.distance < 0 || dst < result.distance))
                        {
                            result.distance = dst;
                            result.u1       = u1;
                            result.u2       = u2;
                        }
                    }
                }

                if (!(dlt < 3) && (result.minDelta < 0 || dlt < result.minDelta))
                    result.minDelta = dlt;
            }

            // --- 2D Newton refinement ---
            {
                double* g = stackalloc double[2];
                double* H = stackalloc double[3];

                short Hsign = NewtonSqDist(&data,
                    ref result.u1, ref result.u2,
                    out result.u1Error, out result.u2Error,
                    out result.distance, out result.distanceError,
                    minrooterr, out result.iterCount2D, 30, g, H);
                Hsign *= (short)(data.a1 * data.a2 >= 0 ? 1 : -1);

                double angerr = EPS * math.PI_DBL * nu;
                result.u1Error += angerr;
                result.u2Error += angerr;

                double tmp     = 2.0 * data.a1 * data.a2;
                double x1      = data.e1 <= 1.0 ? math.cos(result.u1) : math.cosh(result.u1);
                double x2      = data.e2 <= 1.0 ? math.cos(result.u2) : math.cosh(result.u2);
                double tmp1    = data.e1 * x1;
                double tmp2    = data.e2 * x2;
                double r1      = data.a1 * (1.0 - tmp1);
                double r2      = data.a2 * (1.0 - tmp2);
                double rd1     = math.abs(data.a1) * math.sqrt(math.abs(1.0 - tmp1 * tmp1));
                double rd2     = math.abs(data.a2) * math.sqrt(math.abs(1.0 - tmp2 * tmp2));
                double raderr  = math.sqrt(r1 * r1 + r2 * r2) * EPS * nu;
                double graderr = raderr * math.sqrt(rd1 * rd1 + rd2 * rd2) / math.abs(tmp);

                double eigenvalMax = math.abs(H[0] + H[1]) / 2.0
                    + math.sqrt(Sqr((H[0] - H[1]) / 2.0) + H[2] * H[2]);
                double detH  = math.abs(H[0] * H[1] - H[2] * H[2]);
                double uerr3 = eigenvalMax / detH * graderr;
                result.u1Error       += uerr3;
                result.u2Error       += uerr3;
                result.distanceError += eigenvalMax / 2.0 * (Sqr(angerr) + Sqr(graderr) / detH);

                result.distance      *= tmp;
                result.distanceError *= math.abs(tmp);

                result.distance      =  math.sqrt(result.distance);
                result.distanceError += 2.0 * result.distance * raderr + Sqr(raderr);
                result.distanceError /= math.sqrt(Sqr(result.distance) + result.distanceError / 2.0);

                // Diagnostic checks
                result.good = result.good
                    && result.rootCount >= (data.e1 <= 1.0 ? 4 : 2)
                    && result.rootCount % 2 == 0 && result.minDelta > 10
                    && Hsign == +2
                    && math.abs(result.u1Error) <= maxrooterr;
            }

            return result;
        }

        private static double RescaleCoefficients(SAuxData data, Cmplx* c, ref double cerr)
        {
            double znrm = data.e1 <= 1.0
                ? 1.0
                : math.pow(c[0].Norm / c[DEG].Norm, 1.0 / (2 * DEG));
            {
                double tmp = 1.0;
                for (int i = 1; i <= DEG / 2; i++)
                {
                    tmp            *= znrm;
                    c[DEG / 2 + i] *= tmp;
                    c[DEG / 2 - i] /= tmp;
                }

                cerr *= math.max(tmp, 1.0 / tmp);
            }
            return znrm;
        }

        private static void GenerateRootGuess(Cmplx* roots, Cmplx* c, SAuxData data, double znrm, Random rng)
        {
            for (int i = 0; i < DEG; i++)
                roots[i] = 0.0;

            // --- Seed initial root guesses ---
            {
                int    idx        = 4;
                double maxabsroot = MaxRootBound(DEG, c, true) + EPS;
                double minabsroot = data.e1 <= 1.0
                    ? 1.0 / maxabsroot
                    : 1.0 / MaxRootBound(DEG, c, false) + EPS;

                if (data.e1 <= 1.0)
                    roots[0] = Cmplx.Polar(minabsroot, rng.NextDouble() * 2.0 * math.PI_DBL);
                else
                {
                    roots[0] = Cmplx.Polar(1.0 / znrm, 2.0 * math.PI_DBL / 5.0);
                    roots[2] = Cmplx.Polar(1.0 / znrm, data.e2 <= 1.0 ? 3.0 * math.PI_DBL / 5.0 : math.PI_DBL / 2.0);
                }

                // Seed roots near orbital nodes and ±90° from them
                double eta1 = math.sqrt(math.abs(1.0 - data.e1 * data.e1));
                double eta2 = math.sqrt(math.abs(1.0 - data.e2 * data.e2));
                double wPQ1 = math.sqrt(data.P1w * data.P1w + data.Q1w * data.Q1w);
                double wPQ2 = math.sqrt(data.P2w * data.P2w + data.Q2w * data.Q2w);
                double u1__, u2__;

                double* g = stackalloc double[2];
                double* H = stackalloc double[3];

                double u1min = data.e1 <= 0 ? 0 : -math.log(maxabsroot * znrm);
                double u1max = data.e1 <= 0 ? 0 : -math.log(minabsroot * znrm);

                // Ascending node
                double u1asc = Atan2Smart(data.Q1w * eta1, data.e1 * wPQ1 + data.P1w, data.e1 <= 1.0);
                bool   br1   = wPQ1 + data.e1 * data.P1w >= 0;
                double u2asc = Atan2Smart(data.Q2w * eta2, data.e2 * wPQ2 + data.P2w, data.e2 <= 1.0);
                bool   br2   = wPQ2 + data.e2 * data.P2w >= 0;

                // Descending node
                double u1desc = Atan2Smart(-data.Q1w * eta1, data.e1 * wPQ1 - data.P1w, data.e1 <= 1.0);
                bool   br1_   = wPQ1 - data.e1 * data.P1w >= 0;
                double u2desc = Atan2Smart(-data.Q2w * eta2, data.e2 * wPQ2 - data.P2w, data.e2 <= 1.0);
                bool   br2_   = wPQ2 - data.e2 * data.P2w >= 0;

                // Root near ascending node
                {
                    u1__ = u1asc;
                    u2__ = u2asc;
                    SqDist2DIter(&data, br1, br2, ref u1__, ref u2__, out double _, g, H, out double _);
                    if (!(data.e1 <= 1.0)) u1__ = math.clamp(u1__, u1min, u1max);
                    roots[idx++] = Cmplx.Polar(
                        data.e1 <= 1.0 ? 1.0 : math.exp(-u1__) * (br1 ? 1 : -1) / znrm,
                        data.e1 <= 1.0 ? u1__ : rng.NextDouble() * 1e-2 + 1e-3);
                }

                // Root near descending node
                if ((data.e1 <= 1.0 || br1_ != br1) && (data.e2 <= 1.0 || br2_ != br2))
                {
                    u1__ = u1desc;
                    u2__ = u2desc;
                    SqDist2DIter(&data, br1_, br2_, ref u1__, ref u2__, out double _, g, H, out double _);
                    if (!(data.e1 <= 1.0)) u1__ = math.clamp(u1__, u1min, u1max);
                    roots[idx++] = Cmplx.Polar(
                        data.e1 <= 1.0 ? 1.0 : math.exp(-u1__) * (br1_ ? 1 : -1) / znrm,
                        data.e1 <= 1.0 ? u1__ : rng.NextDouble() * 1e-2 + 1e-3);
                }
                else
                {
                    u1__ = data.e1 <= 1.0 || br1_ != br1 ? u1desc : (u1asc + u1desc) / 2.0;
                    u2__ = data.e2 <= 1.0 || br2_ != br2 ? u2desc : (u2asc + u2desc) / 2.0;
                    SqDist2DIter(&data, data.e1 <= 1.0 || !br1, data.e2 <= 1.0 || !br2,
                        ref u1__, ref u2__, out double _, g, H, out double _);
                    if (!(data.e1 <= 1.0)) u1__ = math.clamp(u1__, u1min, u1max);
                    roots[idx++] = Cmplx.Polar(
                        data.e1 <= 1.0 ? 1.0 : math.exp(-u1__) * (br1 ? -1 : 1) / znrm,
                        data.e1 <= 1.0 ? u1__ : rng.NextDouble() * 1e-2 + 1e-3);
                }

                // ±90° from nodes (EE case) or additional node-based seeds (mixed cases)
                if (data.e1 <= 1.0 && data.e2 <= 1.0)
                {
                    // +π/2 from ascending node (in true anomalies)
                    u1__ = math.atan2(-data.P1w * eta1, data.e1 * wPQ1 + data.Q1w);
                    u2__ = math.atan2(-data.P2w * eta2, data.e2 * wPQ2 + data.Q2w);
                    double u2alt = math.atan2(data.P2w * eta2, data.e2 * wPQ2 - data.Q2w);
                    if (2.0 * data.I > math.PI_DBL)
                    {
                        double t = u2__;
                        u2__  = u2alt;
                        u2alt = t;
                    }

                    SqDist2DIter(&data, true, true, ref u1__, ref u2__, out double _, g, H, out double _);
                    roots[idx++] = Cmplx.Polar(1.0, u1__);

                    // -π/2 from ascending node
                    u1__ = math.atan2(data.P1w * eta1, data.e1 * wPQ1 - data.Q1w);
                    SqDist2DIter(&data, true, true, ref u1__, ref u2alt, out double _, g, H, out double _);
                    roots[idx++] = Cmplx.Polar(1.0, u1__);
                }
                else
                {
                    if (data.e1 <= 1.0 || br1_ != br1)
                    {
                        u1__ = u1desc;
                        u2__ = u2asc;
                        SqDist2DIter(&data, br1_, br2, ref u1__, ref u2__, out double _, g, H, out double _);
                        if (!(data.e1 <= 1.0)) u1__ = math.clamp(u1__, u1min, u1max);
                        roots[idx++] = Cmplx.Polar(
                            data.e1 <= 1.0 ? 1.0 : math.exp(-u1__) * (br1_ ? 1 : -1) / znrm,
                            data.e1 <= 1.0 ? u1__ : rng.NextDouble() * 1e-2 + 1e-3);
                    }
                    else
                    {
                        u1__ = (u1asc + u1desc) / 2.0;
                        u2__ = u2asc;
                        SqDist2DIter(&data, !br1, br2, ref u1__, ref u2__, out double _, g, H, out double _);
                        u1__         = math.clamp(u1__, u1min, u1max);
                        roots[idx++] = math.exp(-u1__) * (br1 ? -1 : 1) / znrm;
                    }

                    if (data.e2 <= 1.0 || br2_ != br2)
                    {
                        u1__ = u1asc;
                        u2__ = u2desc;
                        SqDist2DIter(&data, br1, br2_, ref u1__, ref u2__, out double _, g, H, out double _);
                        if (!(data.e1 <= 1.0)) u1__ = math.clamp(u1__, u1min, u1max);
                        roots[idx++] = Cmplx.Polar(
                            data.e1 <= 1.0 ? 1.0 : math.exp(-u1__) * (br1 ? 1 : -1) / znrm,
                            data.e1 <= 1.0 ? u1__ : rng.NextDouble() * 1e-2 + 1e-3);
                    }
                    else
                    {
                        u1__ = u1asc;
                        u2__ = (u2asc + u2desc) / 2.0;
                        SqDist2DIter(&data, br1, !br2, ref u1__, ref u2__, out double _, g, H, out double _);
                        if (!(data.e1 <= 1.0)) u1__ = math.clamp(u1__, u1min, u1max);
                        roots[idx++] = Cmplx.Polar(
                            data.e1 <= 1.0 ? 1.0 : math.exp(-u1__) * (br1 ? 1 : -1) / znrm,
                            data.e1 <= 1.0 ? u1__ : rng.NextDouble() * 1e-2 + 1e-3);
                    }
                }

                // Fill remaining slots
                if (data.e1 <= 1.0)
                    roots[idx] = Cmplx.Polar(minabsroot, rng.NextDouble() * 2.0 * math.PI_DBL);
                else
                {
                    roots[idx]     = Cmplx.Polar(1.0 / znrm, 2.0 * math.PI_DBL / 5.0);
                    roots[idx + 2] = Cmplx.Polar(1.0 / znrm, data.e2 <= 1.0 ? 3.0 * math.PI_DBL / 5.0 : math.PI_DBL / 2.0);
                }
            }
        }
    }
}
