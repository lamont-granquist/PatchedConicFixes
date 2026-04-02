using Unity.Burst;

namespace PatchedConicFixes
{
    /// <summary>
    ///     Result of the MOID computation.
    /// </summary>
    [BurstCompile]
    public struct SMOIDResult
    {
        /// <summary>True if the result is numerically reliable.</summary>
        public bool good;

        /// <summary>Minimum distance between orbits (or -1 if not found).</summary>
        public double distance;

        /// <summary>Numeric uncertainty of the distance.</summary>
        public double distanceError;

        /// <summary>Eccentric anomaly on the first orbit at the MOID point.</summary>
        public double u1;

        /// <summary>Numeric uncertainty of u1.</summary>
        public double u1Error;

        /// <summary>Eccentric anomaly on the second orbit at the MOID point.</summary>
        public double u2;

        /// <summary>Numeric uncertainty of u2.</summary>
        public double u2Error;

        /// <summary>Number of real (or near-real) roots of g(u1) found.</summary>
        public int rootCount;

        /// <summary>
        ///     The minimum delta among all non-real roots.
        ///     Delta measures deviation from the unit circle (elliptic) or real axis (hyperbolic)
        ///     in units of the estimated root error.
        /// </summary>
        public double minDelta;

        /// <summary>Total Newton iterations for 1D polynomial root finding.</summary>
        public ulong iterCount;

        /// <summary>Total Newton iterations for 2D squared-distance refinement.</summary>
        public ulong iterCount2D;

        public static SMOIDResult Default()
        {
            return new SMOIDResult
            {
                good          = true,
                distance      = -1.0,
                distanceError = 0.0,
                u1            = 0.0,
                u1Error       = 0.0,
                u2            = 0.0,
                u2Error       = 0.0,
                rootCount     = 0,
                minDelta      = -1.0,
                iterCount     = 0,
                iterCount2D   = 0
            };
        }
    }
}
