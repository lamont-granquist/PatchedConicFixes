using System.Runtime.CompilerServices;
using Unity.Burst;
using Unity.Mathematics;

namespace PatchedConicFixes
{
    [BurstCompile]
    public static unsafe class PolySolver
    {
        private const double EPS = 2.2204460492503131e-16;

        /// <summary>
        ///     Evaluates polynomial value and first derivative simultaneously via Horner's method.
        ///     Forward: P(z) = c[0] + c[1]z + ... + c[n]z^n  (for |z| ≤ 1)
        ///     Backward: Q(w) = c[n] + c[n-1]w + ... + c[0]w^n  (for |z| > 1, w = 1/z)
        /// </summary>
        public static void PolynomialValDer(int n, Cmplx* c, Cmplx z, out Cmplx val, out Cmplx der, bool forward)
        {
            if (forward)
            {
                val = c[n];
                der = c[n] * n;
                for (int i = n - 1; i > 1; i--)
                {
                    val = val * z + c[i];
                    der = der * z + c[i] * i;
                }

                if (n > 1)
                {
                    val = val * z + c[1];
                    der = der * z + c[1];
                }

                if (n > 0)
                    val = val * z + c[0];
            }
            else
            {
                val = c[0];
                der = c[0] * n;
                for (int i = 1; i < n - 1; i++)
                {
                    val = val * z + c[i];
                    der = der * z + c[i] * (n - i);
                }

                if (n > 1)
                {
                    val = val * z + c[n - 1];
                    der = der * z + c[n - 1];
                }

                if (n > 0)
                    val = val * z + c[n];
            }
        }

        /// <summary>
        /// Evaluates polynomial value, first derivative, and second derivative simultaneously
        /// via Horner's method.
        ///
        /// Forward (|z| ≤ 1): P(z) = c[0] + c[1]z + ... + c[n]z^n
        /// Backward (|z| > 1): Q(w) = c[n] + c[n-1]w + ... + c[0]w^n, where w = 1/z
        ///
        /// Used by RootError for the quadratic approximation of root uncertainty.
        /// </summary>
        public static void PolynomialValDer2(int n, Cmplx* c, Cmplx z,
            out Cmplx val, out Cmplx der, out Cmplx der2,
            bool forward)
        {
            if (forward)
            {
                val  = c[n];
                der  = c[n] * n;
                der2 = c[n] * (double)(n * (n - 1));
                for (int i = n - 1; i > 2; i--)
                {
                    val  = val  * z + c[i];
                    der  = der  * z + c[i] * i;
                    der2 = der2 * z + c[i] * (double)(i * (i - 1));
                }
                if (n > 2)
                {
                    val  = val  * z + c[2];
                    der  = der  * z + c[2] * 2.0;
                    der2 = der2 * z + c[2] * 2.0;
                }
                if (n > 1)
                {
                    val = val * z + c[1];
                    der = der * z + c[1];
                }
                if (n > 0)
                    val = val * z + c[0];
            }
            else
            {
                val  = c[0];
                der  = c[0] * n;
                der2 = c[0] * (double)(n * (n - 1));
                for (int i = 1; i < n - 2; i++)
                {
                    val  = val  * z + c[i];
                    der  = der  * z + c[i] * (double)(n - i);
                    der2 = der2 * z + c[i] * (double)((n - i) * (n - i - 1));
                }
                if (n > 2)
                {
                    val  = val  * z + c[n - 2];
                    der  = der  * z + c[n - 2] * 2.0;
                    der2 = der2 * z + c[n - 2] * 2.0;
                }
                if (n > 1)
                {
                    val = val * z + c[n - 1];
                    der = der * z + c[n - 1];
                }
                if (n > 0)
                    val = val * z + c[n];
            }
        }

        /// <summary>
        ///     Computes f(z)/f'(z) for a polynomial of degree n.
        ///     Returns zero if the numerator is zero (exact root).
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx Ratio1(int n, Cmplx* c, Cmplx z, bool forward)
        {
            PolynomialValDer(n, c, z, out Cmplx f, out Cmplx fd, forward);
            return f == 0.0 ? f : f / fd;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Cmplx Ratio2(Cmplx* c, Cmplx z, bool forward)
        {
            Cmplx num, den;

            if (forward)
            {
                num = (c[2] * z + c[1]) * z + c[0];
                den = 2.0 * c[2] * z + c[1];
            }
            else
            {
                num = (c[0] * z + c[1]) * z + c[2];
                den = 2.0 * c[0] * z + c[1];
            }

            return num == 0.0 ? 0.0 : num / den;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static Cmplx Ratio4(Cmplx* c, Cmplx z, bool forward)
        {
            Cmplx num, den;

            if (forward)
            {
                num = (((c[4] * z + c[3]) * z + c[2]) * z + c[1]) * z + c[0];
                den = ((4.0 * c[4] * z + 3.0 * c[3]) * z + 2.0 * c[2]) * z + c[1];
            }
            else
            {
                num = (((c[0] * z + c[1]) * z + c[2]) * z + c[3]) * z + c[4];
                den = ((4.0 * c[0] * z + 3.0 * c[1]) * z + 2.0 * c[2]) * z + c[3];
            }

            return num == 0.0 ? 0.0 : num / den;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void PolishRoot(Cmplx* c, ref Cmplx x)
        {
            if (x.Norm <= 1.0)
                x -= Ratio4(c, x, true);
            else
                x /= 1.0 - x * Ratio4(c, x.Inverse(), false);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static void PolishRoot2(Cmplx* c, ref Cmplx x)
        {
            if (x.Norm <= 1.0)
                x -= Ratio2(c, x, true);
            else
                x /= 1.0 - x * Ratio2(c, x.Inverse(), false);
        }

        /// <summary>
        ///     Solves c[0] + c[1]z + c[2]z² = 0.
        ///     Writes two roots into roots[0] and roots[1].
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static int SolveQuadratic(Cmplx* c, Cmplx* roots)
        {
            var D = Cmplx.Sqrt(c[1] * c[1] - 4.0 * c[2] * c[0]);

            double dot = c[1].Re * D.Re + c[1].Im * D.Im;
            Cmplx  tmp = dot >= 0 ? -c[1] - D : -c[1] + D;

            roots[0] = tmp / (2.0 * c[2]);
            roots[1] = 2.0 * c[0] / tmp;

            PolishRoot2(c, ref roots[0]);
            PolishRoot2(c, ref roots[1]);

            return 2;
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        private static double Cbrt(double x)
        {
            // FIXME: this could be optimized even more for speed + accuracy with one of the magic constant cbrt() methods
            if (x == 0.0) return x;

            double y = math.sign(x) * math.exp(math.log(math.abs(x)) / 3.0);
            double t = y * y * y;
            return y - y * (t - x) / (2.0 * t + x);
        }

        /// <summary>
        ///     Solves c[0] + c[1]z + c[2]z² + c[3]z³ + c[4]z⁴ = 0.
        ///     Writes four roots into roots[0..3].
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static int SolveQuartic(Cmplx* c, Cmplx* roots)
        {
            Cmplx r = 4.0 * c[4];
            Cmplx p = 2.0 * r * c[2] - 3.0 * Cmplx.Sqr(c[3]);
            Cmplx q = 2.0 * c[3] * (Cmplx.Sqr(c[3]) - r * c[2]) + Cmplx.Sqr(r) * c[1];

            Cmplx D0 = Cmplx.Sqr(c[2]) - 3.0 * c[3] * c[1] + 3.0 * c[0] * r;
            Cmplx D1 = 2.0 * c[2] * Cmplx.Sqr(c[2])
                - 9.0 * c[3] * c[2] * c[1]
                + 27.0 * Cmplx.Sqr(c[3]) * c[0]
                + 27.0 * Cmplx.Sqr(c[1]) * c[4]
                - 18.0 * r * c[2] * c[0];

            var    tmp = Cmplx.Sqrt(Cmplx.Sqr(D1) - 4.0 * D0 * Cmplx.Sqr(D0));
            double dot = D1.Re * tmp.Re + D1.Im * tmp.Im;
            Cmplx  Qc  = (dot >= 0 ? D1 + tmp : D1 - tmp) / 2.0;

            double Qabs = Cbrt(Qc.Abs);
            double Qarg = Qc.Phase / 3.0;

            const double TWO_PI_OVER3 = 2.0 * math.PI_DBL / 3.0;

            var Q1 = Cmplx.Polar(Qabs, Qarg);
            var Q2 = Cmplx.Polar(Qabs, Qarg + TWO_PI_OVER3);
            var Q3 = Cmplx.Polar(Qabs, Qarg - TWO_PI_OVER3);

            Cmplx D0_over_Q1 = D0 == 0.0 ? 0.0 : D0 / Q1;
            Cmplx D0_over_Q2 = D0 == 0.0 ? 0.0 : D0 / Q2;
            Cmplx D0_over_Q3 = D0 == 0.0 ? 0.0 : D0 / Q3;

            Cmplx val1 = (Q1 + D0_over_Q1) * r - p;
            Cmplx val2 = (Q2 + D0_over_Q2) * r - p;
            Cmplx val3 = (Q3 + D0_over_Q3) * r - p;

            double n1 = val1.Norm, n2 = val2.Norm, n3 = val3.Norm;

            Cmplx bestVal;
            if (n1 > n2 && n1 > n3) bestVal = val1;
            else if (n2 > n3) bestVal       = val2;
            else bestVal                    = val3;

            tmp = Cmplx.Sqrt(bestVal / 3.0);
            Cmplx S = c[4].Re >= 0.0 ? tmp : -tmp;

            Cmplx q_over_S = q == 0.0 ? 0.0 : q / S;

            tmp      = Cmplx.Sqrt(q_over_S - Cmplx.Sqr(S) - p);
            roots[0] = (-c[3] - S + tmp) / r;
            roots[1] = (-c[3] - S - tmp) / r;

            tmp      = Cmplx.Sqrt(-q_over_S - Cmplx.Sqr(S) - p);
            roots[2] = (-c[3] + S + tmp) / r;
            roots[3] = (-c[3] + S - tmp) / r;

            PolishRoot(c, ref roots[0]);
            PolishRoot(c, ref roots[1]);
            PolishRoot(c, ref roots[2]);
            PolishRoot(c, ref roots[3]);

            return 4;
        }

        /// <summary>
        ///     Divides out a known root from the polynomial using Ruffini-Horner deflation.
        ///     For numerical stability, extracts (z - root) when |root| ≤ 1,
        ///     or (1/z - 1/root) when |root| > 1.
        ///     Returns a pointer to the first coefficient of the reduced polynomial.
        /// </summary>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx* ExtractLinearFactor(int n, Cmplx* c, Cmplx root)
        {
            double rootNorm = root.Norm;
            if (rootNorm <= 1.0)
            {
                for (int i = n - 1; i > 0; i--)
                    c[i] += c[i + 1] * root;
                return c + 1;
            }

            Cmplx iroot = root.Conjugate() / rootNorm;
            for (int i = 1; i < n; i++)
                c[i] += c[i - 1] * iroot;
            return c;
        }

        /// <summary>
        ///     Finds one complex root of polynomial c[0] + c[1]z + ... + c[n]z^n via Newton's method.
        ///     Two-phase convergence:
        ///     Phase 1: Converge to maxeps relative tolerance, with cycle detection and random restarts.
        ///     Phase 2: Push toward mineps or until no further progress (machine precision).
        ///     Inner/outer switching: evaluates in z when |z| is small, in 1/z when |z| is large,
        ///     for numerical stability.
        ///     Returns total iteration count.
        /// </summary>
        public static ulong Newton(int n, Cmplx* c, ref Cmplx z, double maxeps, double mineps, ulong maxcount, ulong restartcount, Random rng)
        {
            maxeps = math.abs(maxeps);
            mineps = math.abs(mineps);

            // Enforce minimum tolerance floors
            if (!(maxeps >= 2.0 * mineps)) maxeps = 2.0 * mineps;
            if (!(maxeps >= 2.0 * EPS)) maxeps    = 2.0 * EPS;

            Cmplx z1 = z, z2 = z1;
            Cmplx dz;

            ulong  iterCnt = 0;
            double err2    = 0.0, nrm;
            bool   inner   = true;

            // --- Phase 1: converge to maxeps with cycle detection ---
            do
            {
                iterCnt++;

                // Detect cycling: z ≈ z2 (period-2 orbit) or hit iteration budget
                bool cycled = (iterCnt > 2 && (z - z2).Norm < (z - z1).Norm * 1e-6)
                    || iterCnt % maxcount == 0;

                Cmplx randMlt = default;
                if (cycled)
                {
                    double rMag = (2.0 * rng.NextDouble() + 1.0) / 3.0;
                    double rArg = rng.NextDouble() * 2.0 * math.PI_DBL;
                    randMlt = Cmplx.Polar(rMag, rArg);
                }

                nrm = z.Norm;

                // Hysteresis on inner/outer switching to avoid thrashing at |z| ≈ 1
                if (inner)
                {
                    if (nrm > 5.0) inner = false;
                }
                else
                {
                    if (2.0 * nrm < 1.0) inner = true;
                }

                z2 = z1;
                z1 = z;

                if (inner)
                {
                    dz =  Ratio1(n, c, z, true);
                    z  -= cycled ? dz * randMlt : dz;
                }
                else
                {
                    nrm =  1.0 / nrm;
                    dz  =  Ratio1(n, c, z.Conjugate() * nrm, false); // dw, not dz
                    z   /= 1.0 - z * (cycled ? dz * randMlt : dz);
                }

                err2 = dz.Norm;
            } while (err2 > maxeps * maxeps * nrm && iterCnt < maxcount * restartcount);

            // --- Phase 2: push toward mineps or machine precision ---
            if (err2 > mineps * mineps * nrm)
            {
                ulong  iterCnt2 = 0;
                double derr2;
                do
                {
                    nrm = z.Norm;
                    if (inner)
                    {
                        dz =  Ratio1(n, c, z, true);
                        z  -= dz;
                    }
                    else
                    {
                        nrm =  1.0 / nrm;
                        dz  =  Ratio1(n, c, z.Conjugate() * nrm, false);
                        z   /= 1.0 - z * dz;
                    }

                    double err2New = dz.Norm;
                    derr2 = err2 - err2New;
                    err2  = err2New;
                    iterCnt2++;
                } while (derr2 > EPS * err2
                         && err2 > mineps * mineps * nrm
                         && iterCnt2 <= maxcount);

                iterCnt += iterCnt2;
            }

            return iterCnt;
        }

        /// <summary>
        ///     Finds all complex roots of a degree-n polynomial.
        ///     P(z) = c[0] + c[1]z + ... + c[n]z^n = 0
        ///     Uses Newton's method with deflation for roots n down to 5,
        ///     then finishes with a closed-form quartic solver.
        ///     The roots array must be pre-allocated with at least n elements.
        ///     Non-zero entries in roots[] are used as initial guesses for Newton;
        ///     zero entries get a predicted starting guess based on the previous root
        ///     and the polynomial's symmetry structure.
        ///     For trigonometric polynomials (from the MOID problem): if z is a root,
        ///     then 1/z* is also a root, so the next guess is set to z/|z|² = 1/z*.
        ///     For real-coefficient polynomials: if z is a root, then z* is also a root,
        ///     so the next guess is set to conj(z).
        ///     Returns total Newton iteration count.
        /// </summary>
        public static ulong PolynomialRoots(int n, Cmplx* c, Cmplx* roots, double minerr, double maxerr, bool trigonometric, Random rng)
        {
            Cmplx* c_      = c;
            ulong  iterCnt = 0;

            // Extract roots one at a time until degree 4 remains
            for (int i = 0, degree = n; degree > 4; i++)
            {
                iterCnt += Newton(degree, c_, ref roots[i], maxerr, minerr, 300, 10, rng);
                c_      =  ExtractLinearFactor(degree, c_, roots[i]);
                degree--;

                // Seed the next root's initial guess if it hasn't been pre-set
                if (degree > 4 && roots[i + 1] == 0.0)
                {
                    if (trigonometric)
                    {
                        // Trig polynomial symmetry: z is a root => 1/z* is a root
                        roots[i + 1] = roots[i] / roots[i].Norm;
                    }
                    else
                    {
                        // Real-coefficient symmetry: z is a root => z* is a root
                        if (math.abs(roots[i].Im) > math.abs(roots[i].Re) * 1e-3)
                            roots[i + 1] = roots[i].Conjugate();
                        else
                            roots[i + 1] = new Cmplx(roots[i].Im, roots[i].Re);
                    }
                }
            }

            // Solve the remaining quartic
            iterCnt += (ulong)SolveQuartic(c_, roots + n - 4);

            return iterCnt;
        }

        // XXX: for tests
        public static Random MakeRng() => new Random(12345);
    }
}
