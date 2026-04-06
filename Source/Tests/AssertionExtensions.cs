using System;
using System.Globalization;
using Xunit.Sdk;
using static PatchedConicFixes.Statics;
using static System.Math;

namespace PatchedConicFixes.Tests
{
    public static class AssertionExtensions
    {
        public static void ShouldEqual(this int actual, int expected)
        {
            if (actual != expected)
                throw new XunitException(
                    string.Format(CultureInfo.CurrentCulture, "Expected integer to be '{0}', but was '{1}'", expected, actual)
                );
        }

        // A rtol==atol comparison between double precision floats
        public static void ShouldEqual(this double actual, double expected, double epsilon = EPS)
        {
            if (double.IsNaN(epsilon) || double.IsNegativeInfinity(epsilon) || epsilon < 0.0)
                throw new ArgumentException("Epsilon must be greater than or equal to zero", nameof(epsilon));

            if (!NearlyEqual(actual, expected, epsilon))
                throw new ApproximateEqualException(
                    string.Format(CultureInfo.CurrentCulture, "{0:G17}", expected),
                    string.Format(CultureInfo.CurrentCulture, "{0:G17}", actual),
                    epsilon
                );
        }
    }
}
