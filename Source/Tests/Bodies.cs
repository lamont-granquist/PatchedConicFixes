using System;
using System.Collections.Generic;
using System.Runtime.Serialization;

namespace PatchedConicFixes.Tests
{
    public static class Bodies
    {
        public static (CelestialBody earth, CelestialBody moon) MakeEarthMoon()
        {
            return MakeParentChild(3.9860043543609598e+14, 924649202.461023, 4.9028000661637961e+12, 66167158.6569544, 28.3626779079849, 0.0532814935368257, 384308437.770707, 2.29661616112602, 199.764093016082, 3.88686980063246, -31542641.784);
        }

        public static (CelestialBody parent, CelestialBody child) MakeParentChild(
            double parentMu, double parentSoi,
            double childMu, double childSoi,
            double inc, double e, double sma, double lan, double argPe, double mEp, double epoch)
        {
            var parent = MakeBody(parentMu, parentSoi);
            var child = MakeBody(childMu, childSoi);

            if (!(FormatterServices.GetUninitializedObject(typeof(OrbitDriver)) is OrbitDriver orbitDriver))
                throw new InvalidOperationException("Failed to create OrbitDriver");
            orbitDriver.orbit         = new Orbit(inc, e, sma, lan, argPe, mEp, epoch, parent);
            orbitDriver.celestialBody = child;
            child.orbitDriver         = orbitDriver;

            parent.orbitingBodies.Add(child);

            return (parent, child);
        }

        private static CelestialBody MakeBody(double mu, double soi)
        {
            if (!(FormatterServices.GetUninitializedObject(typeof(CelestialBody)) is CelestialBody cb))
                throw new InvalidOperationException("Failed to create CelestialBody");

            cb.gravParameter = mu;
            cb.sphereOfInfluence = soi;
            cb.orbitingBodies = new List<CelestialBody>();

            return cb;
        }
    }
}
