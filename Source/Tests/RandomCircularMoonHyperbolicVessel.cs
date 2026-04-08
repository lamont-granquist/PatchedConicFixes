using System;
using System.Collections.Generic;
using Xunit;
using Xunit.Abstractions;

namespace PatchedConicFixes.Tests
{
    [Collection("KSP")]
    public class RandomCircularMoonHyperbolicVessel
    {
        private readonly ITestOutputHelper _output;

        public RandomCircularMoonHyperbolicVessel(ITestOutputHelper output)
        {
            _output = output;
        }

        // 15790, 12279, 425 hangs stock?
        public static IEnumerable<object[]> Seeds()
        {
            for (int i = 0; i <= 500; i++)
                yield return new object[] { i };
        }

        /// <summary>
        ///     Constructs a guaranteed encounter by working backwards:
        ///     1. Create a moon on a circular equatorial orbit around a parent body
        ///     2. Pick a random point on the moon's SOI sphere at a known encounter time
        ///     3. Construct a vessel orbit passing through that SOI point with an infalling velocity
        ///     4. Rewind time and run the encounter solver
        ///     5. Assert the solver finds the encounter (or an earlier valid one)
        ///     For each geometry, 10 different rewind amounts are tested to exercise different
        ///     solver starting configurations.
        /// </summary>
        [Theory]
        [MemberData(nameof(Seeds))]
        public void RandomlyConstructedEncounter(int seed)
        {
            Logger.Register(o => _output.WriteLine((string)o));
            var rng = new Random(seed);

            // build a random system
            double parentMu   = LogUniform(rng, 1e10, 1e15);
            double moonSma    = LogUniform(rng, 1e7, 5e8);
            double moonMu     = LogUniform(rng, 1e8, 1e12);
            double foo        = Uniform(rng, 0.01, 0.08);
            double moonSoi    = moonSma * foo;
            double moonMEpoch = Uniform(rng, 0, 2.0 * Math.PI);

            (CelestialBody parent, CelestialBody moon) = Bodies.MakeParentChild(
                parentMu, 10 * moonSma,
                moonMu, moonSoi,
                0, 0, moonSma, 0, 0, moonMEpoch, 0
            );

            /*
            Logger.Print("Bodies.MakeParentChild(");
            Logger.Print($"  {parentMu}, {10 * moonSma}");
            Logger.Print($"  {moonMu}, {moonSoi}");
            Logger.Print($"  0, 0, {moonSma}, 0, 0, {moonMEpoch}, 0");
            Logger.Print(")");
            */

            // Pick the SOI encounter time
            double moonPeriod = 2.0 * Math.PI * Math.Sqrt(moonSma * moonSma * moonSma / parentMu);
            double tEnc       = moonPeriod * Uniform(rng, 10.0, 50.0);

            Logger.Print($"ENCOUNTER TIME: {tEnc}");

            // Moon state at encounter time
            Orbit    moonOrbit = moon.orbitDriver.orbit;
            Vector3d moonPos   = moonOrbit.getRelativePositionAtUT(tEnc);
            Vector3d moonVel   = moonOrbit.getOrbitalVelocityAtUT(tEnc);

            // Vessel position at encounter time
            Vector3d soiOffset = moonSoi * RandomUnitVector(rng);
            Vector3d vesselPos = moonPos + soiOffset;
            double   vesselR   = vesselPos.magnitude;

            // Vessel velocity at encounter time
            double vCirc = Math.Sqrt(parentMu / vesselR);

            Vector3d vesselVel = Vector3d.zero;

            for (int attempt = 0; attempt < 200; attempt++)
            {
                double   speed        = vCirc * Uniform(rng, 1.42, 3.0);
                Vector3d candidateVel = speed * RandomUnitVector(rng);

                if (Vector3d.Dot(soiOffset, candidateVel - moonVel) < 0)
                {
                    vesselVel = candidateVel;
                    break;
                }
            }

            Assert.False(vesselVel == Vector3d.zero, "Failed to generate infalling velocity after 200 attempts");

            // Construct Orbit from state vectors at encounter
            var vesselOrbit = new Orbit();
            vesselOrbit.UpdateFromStateVectors(vesselPos, vesselVel, parent, tEnc);
            Logger.Print($"{vesselOrbit.referenceBody.gravParameter}");
            Logger.Print($"{vesselOrbit.meanMotion}");
            double nuStart = -vesselOrbit.TrueAnomalyAtRadius(10 * moonSma);
            double nuEnd = vesselOrbit.TrueAnomalyAtUT(tEnc);

            /*
            _output.WriteLine($"Seed {seed}: parentMu={parentMu:E3} moonSma={moonSma:E3} moonSoi={moonSoi:E3}");
            _output.WriteLine($"  vesselOrbit: e={vesselOrbit.eccentricity:F6} sma={vesselOrbit.semiMajorAxis:E3} period={vesselPeriod:F1}");
            _output.WriteLine($"  tEnc={tEnc:F3} moonPeriod={moonPeriod:F1}");


            Logger.Print($"{Helpers.OrbitDataString(moonOrbit)}");
            Logger.Print($"{Helpers.OrbitDataString(vesselOrbit)}");
            */

            // --- Try 10 different rewind amounts ---
            int failures = 0;

            for (int trial = 0; trial < 10; trial++)
            {
                // Rewind by 10%-90% of the vessel's period
                double rewindFraction = Uniform(rng, 0.1, 0.9);
                double nuRewind       = nuStart + (nuEnd - nuStart) * rewindFraction;
                double rewind = tEnc + vesselOrbit.GetDTforTrueAnomaly(nuRewind, 0);
                double startEpoch     = tEnc - rewind;

                vesselOrbit.GetOrbitalStateVectorsAtUT(startEpoch, out Vector3d pos, out Vector3d vel);

                /*
                // Build a fresh orbit for each trial (CheckEncounter mutates it)
                Logger.Print($"startEpoch = {startEpoch}");
                Logger.Print($"p.UpdateFromStateVectors(new Vector3d({pos.x}, {pos.y}, {pos.z}), new Vector3d({vel.x}, {vel.y}, {vel.z}), parent, startEpoch)");
                Logger.Print($"p.StartUT = {startEpoch}");
                Logger.Print($"p.EndUT = {startEpoch + vesselPeriod}");
                */

                var p = new Orbit();
                p.UpdateFromStateVectors(pos, vel, parent, startEpoch);
                p.StartUT = startEpoch;
                p.EndUT   = tEnc + vesselOrbit.GetDTforTrueAnomaly(-nuStart, 0);;

                var nextPatch = new Orbit();
                var pars      = new PatchedConics.SolverParameters();

                //bool result = PatchedConics.CheckEncounter(p, nextPatch, startEpoch, moon.orbitDriver, moon, pars, false); // test against actual stock
                //bool result = HarmonyPatches.CheckEncounterOldMOID(p, nextPatch, startEpoch, moon.orbitDriver, moon, pars, false); // test against old MOID code
                bool result = HarmonyPatches.CheckEncounter(p, nextPatch, startEpoch, moon.orbitDriver, moon, pars, false);

                if (!result)
                {
                    _output.WriteLine($"  Trial {trial}: FAIL - no encounter found (rewind={rewind:F1}, startEpoch={startEpoch:F3})");
                    failures++;
                    continue;
                }

                // Encounter found — validate it
                Vector3d vesselPosAtSOI = p.getPositionAtUT(p.EndUT);
                Vector3d moonPosAtSOI   = moonOrbit.getPositionAtUT(p.EndUT);
                double   distAtSOI      = (vesselPosAtSOI - moonPosAtSOI).magnitude;

                // The SOI crossing distance should match the moon's SOI
                bool soiDistOk = Math.Abs(distAtSOI - moonSoi) / moonSoi < 0.01; // 1% tolerance

                if (!soiDistOk)
                {
                    _output.WriteLine($"  Trial {trial}: FAIL - SOI distance mismatch (dist={distAtSOI:F1}, soi={moonSoi:F1})");
                    failures++;
                    continue;
                }

                // Check timing: did we find our encounter or an earlier one?
                // NB: if we find an EARLIER encounter by rewinding that is considered good enough
                if (p.EndUT > tEnc + 11.0)
                {
                    // We should never find a FIRST encounter LATER than the one we've constructed to happen.
                    _output.WriteLine($"  Trial {trial}: FAIL - found later encounter (found={p.EndUT:F3}, expected<={tEnc:F3})");
                    failures++;
                    continue;
                }

                string timing = Math.Abs(p.EndUT - tEnc) < 10.0 ? "exact" : "earlier";
                _output.WriteLine($"  Trial {trial}: OK ({timing}, endUT={p.EndUT:F3}, tEnc={tEnc:F3})");
            }

            Assert.True(failures == 0, $"Seed {seed}: {failures}/10 trials failed");
        }

        #region Helpers

        private static double Uniform(Random rng, double lo, double hi) => lo + (hi - lo) * rng.NextDouble();

        private static double LogUniform(Random rng, double lo, double hi)
        {
            double logLo = Math.Log(lo);
            double logHi = Math.Log(hi);
            return Math.Exp(logLo + (logHi - logLo) * rng.NextDouble());
        }

        private static Vector3d RandomUnitVector(Random rng)
        {
            double z   = 2.0 * rng.NextDouble() - 1.0;
            double phi = 2.0 * Math.PI * rng.NextDouble();
            double r   = Math.Sqrt(1.0 - z * z);
            return new Vector3d(r * Math.Cos(phi), r * Math.Sin(phi), z);
        }

        #endregion
    }
}
