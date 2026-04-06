using System.Runtime.CompilerServices;
using static System.Math;

namespace PatchedConicFixes
{
    public static class Statics
    {
        public const double TAU = 2 * PI;

        /// <summary>
        ///     Normal machine epsilon.  The Double.Epsilon in C# is one ULP above zero which is somewhat useless.
        /// </summary>
        public const double EPS = 2.2204460492503131e-16;

        /// <summary>
        ///     Twice machine epsilon.
        /// </summary>
        public const double EPS2 = EPS * 2;

        public const double DEG2RAD = PI / 180.0;
        public const double RAD2DEG = 180.0 / PI;

        /// <summary>
        ///     Convert Degrees to Radians.
        /// </summary>
        /// <param name="deg">degrees</param>
        /// <returns>radians</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Deg2Rad(double deg) => deg * DEG2RAD;

        /// <summary>
        ///     Convert Radians to Degrees.
        /// </summary>
        /// <param name="rad">Radians</param>
        /// <returns>Degrees</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static double Rad2Deg(double rad) => rad * RAD2DEG;

        /// <summary>
        ///     Helper to check if a value is finite (not NaN or Ininity).
        /// </summary>
        /// <param name="x">Value</param>
        /// <returns>True if the value is finite</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool IsFinite(double x) => !double.IsNaN(x) && !double.IsInfinity(x);

        /// <summary>
        ///     Compares two double values with relative and absolute tolerance.
        /// </summary>
        /// <param name="num">first value</param>
        /// <param name="reference">reference value</param>
        /// <param name="rtol">relative tolerance (e.g. 1e-15)</param>
        /// <param name="atol">absolute tolerance (e.g. 1e-15)</param>
        /// <param name="equalNan">if set, treats two NaN values as equal</param>
        /// <returns>true if the values are nearly the same</returns>
        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool NearlyEqual(double num, double reference, double rtol, double atol, bool equalNan = false)
        {
            if (num.Equals(reference))
                return true;

            if (equalNan && double.IsNaN(num) && double.IsNaN(reference))
                return true;

            if (!IsFinite(num) || !IsFinite(reference))
                return false;

            // see scipy.isclose() - symmetric in arguments
            return Abs(num - reference) <= atol + rtol * Max(Abs(num), Abs(reference));
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool NearlyEqual(double num, double reference, double tol = EPS) => NearlyEqual(num, reference, tol, tol);
    }
}
