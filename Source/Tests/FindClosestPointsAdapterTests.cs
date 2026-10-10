using System;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using HarmonyLib;
using Xunit;
using Xunit.Abstractions;
using Random = System.Random;

namespace PatchedConicFixes.Tests
{
    [Collection("KSP")]
    public class FindClosestPointsAdapterTests
    {
        private readonly ITestOutputHelper _output;

        public FindClosestPointsAdapterTests(ITestOutputHelper output)
        {
            _output = output;
            Logger.Register(o => _output.WriteLine((string)o));
        }

        private struct Result
        {
            public int    Count;
            public double CD, CCD, FFp, FFs, SFp, SFs;
        }

        private static Result Call(Orbit.FindClosestPointsDelegate solver, Orbit p, Orbit s)
        {
            var r          = new Result();
            int iterations = 0;
            r.Count = solver(p, s, ref r.CD, ref r.CCD, ref r.FFp, ref r.FFs, ref r.SFp, ref r.SFs, 0.0001, 20, ref iterations);
            return r;
        }

        // The two crossing points of a coplanar ellipse and circle with the same semi-major axis are where
        // r = a(1-e^2)/(1+e*cos(v)) = a, so cos(v) = -e.
        [Fact]
        public void CoplanarEllipseAndCircleReturnTheAnalyticCrossings()
        {
            (CelestialBody earth, _) = Bodies.MakeEarthMoon();

            const double sma = 7000000;
            const double ecc = 0.2;

            var ellipse = new Orbit(0, ecc, sma, 0, 0, 0, 0, earth);
            var circle  = new Orbit(0, 0, sma, 0, 0, 0, 0, earth);

            Result r = Call(FindClosestPointsAdapter.FindClosestPoints, ellipse, circle);

            Assert.Equal(2, r.Count);
            Assert.Equal(Math.Acos(-ecc), Math.Abs(r.FFp), 7);
            Assert.Equal(Math.Acos(-ecc), Math.Abs(r.SFp), 7);
            Assert.True(r.FFp * r.SFp < 0);
            Assert.InRange(r.CD, 0, 1e-3);
            Assert.InRange(r.CCD, 0, 1e-3);
        }

        [Fact]
        public unsafe void OneMinimumIsCopiedIntoBothSlots()
        {
            Baluev.MoidInfo* info = stackalloc Baluev.MoidInfo[1];
            info[0] = new Baluev.MoidInfo { dst = 1, ta1 = 2, ta2 = 3 };

            double CD = 0, CCD = 0, FFp = 0, FFs = 0, SFp = 0, SFs = 0;
            int    count = FindClosestPointsAdapter.CopyMinima(info, 0, -1, ref CD, ref CCD, ref FFp, ref FFs, ref SFp, ref SFs);

            Assert.Equal(1, count);
            Assert.Equal((1.0, 2.0, 3.0), (CD, FFp, FFs));
            Assert.Equal((1.0, 2.0, 3.0), (CCD, SFp, SFs));
        }

        [Fact]
        public unsafe void NoMinimaLeavesTheOutputsUntouched()
        {
            Baluev.MoidInfo* info = stackalloc Baluev.MoidInfo[1];

            double CD = 1, CCD = 2, FFp = 3, FFs = 4, SFp = 5, SFs = 6;
            int    count = FindClosestPointsAdapter.CopyMinima(info, -1, -1, ref CD, ref CCD, ref FFp, ref FFs, ref SFp, ref SFs);

            Assert.Equal(0, count);
            Assert.Equal((1.0, 2.0, 3.0, 4.0, 5.0, 6.0), (CD, CCD, FFp, FFs, SFp, SFs));
        }

        private static IEnumerable<(Orbit p, Orbit s)> RandomOrbitPairs(CelestialBody body, int seed, int n)
        {
            var rng = new Random(seed);

            Orbit RandomOrbit() => new Orbit(180 * rng.NextDouble(), 0.9 * rng.NextDouble(), 7000000 + 50000000 * rng.NextDouble(),
                360 * rng.NextDouble(), 360 * rng.NextDouble(), 2 * Math.PI * rng.NextDouble(), 0, body);

            for (int i = 0; i < n; i++)
                yield return (RandomOrbit(), RandomOrbit());
        }

        // The stock solver returns two of the minima it finds, in no particular order, so the best minimum returned by the adapter
        // (the global MOID among same-side minima) should never be worse than either of them.
        [Fact]
        public void NeverWorseThanTheStockSolver()
        {
            (CelestialBody earth, _) = Bodies.MakeEarthMoon();

            int compared = 0;

            foreach ((Orbit p, Orbit s) in RandomOrbitPairs(earth, 42, 500))
            {
                Result stock   = Call(Orbit._FindClosestPoints, p, s);
                Result adapter = Call(FindClosestPointsAdapter.FindClosestPoints, p, s);

                if (stock.Count < 1)
                    continue;

                Assert.True(adapter.Count >= 1, $"no minima for {HarmonyPatches.OrbitString(p)} / {HarmonyPatches.OrbitString(s)}");

                double stockBest = stock.Count > 1 ? Math.Min(stock.CD, stock.CCD) : stock.CD;

                Assert.True(adapter.CD <= stockBest * (1 + 1e-6) + 1e-3,
                    $"adapter {adapter.CD} > stock {stockBest} for {HarmonyPatches.OrbitString(p)} / {HarmonyPatches.OrbitString(s)}");

                if (adapter.Count > 1)
                    Assert.True(adapter.CD <= adapter.CCD);

                compared++;
            }

            _output.WriteLine($"compared {compared} orbit pairs");
            Assert.True(compared > 400);
        }

        [Fact]
        public void InstallOnlyReplacesTheStockDelegate()
        {
            Orbit.FindClosestPointsDelegate saved = Orbit.FindClosestPoints;

            try
            {
                Orbit.FindClosestPoints = Orbit._FindClosestPoints;
                FindClosestPointsAdapter.Install();
                Assert.Equal(typeof(FindClosestPointsAdapter), Orbit.FindClosestPoints.Method.DeclaringType);

                // installing twice is harmless
                FindClosestPointsAdapter.Install();
                Assert.Equal(typeof(FindClosestPointsAdapter), Orbit.FindClosestPoints.Method.DeclaringType);

                // another mod's replacement is left alone
                Orbit.FindClosestPointsDelegate other = InterceptDisplay.SelectDisplayIntercepts;
                Orbit.FindClosestPoints = other;
                FindClosestPointsAdapter.Install();
                Assert.Same(other, Orbit.FindClosestPoints);
            }
            finally
            {
                Orbit.FindClosestPoints = saved;
            }
        }

        [Fact]
        public void DisplayInterceptsAcrossDifferentBodiesAreEmpty()
        {
            (CelestialBody earth, CelestialBody moon) = Bodies.MakeEarthMoon();

            var p = new Orbit(0, 0.2, 7000000, 0, 0, 0, 0, earth);
            var s = new Orbit(0, 0, 7000000, 0, 0, 0, 0, moon);

            Result r = Call(InterceptDisplay.SelectDisplayIntercepts, p, s);

            Assert.Equal(0, r.Count);
            Assert.Equal(0.0, r.FFp);
            Assert.Equal(0.0, r.SFp);
        }

        // The display intercepts are the two minima with the smallest separation at the time that p passes them, ordered by
        // separation, and that time is computed exactly as the stock callers compute it.
        [Fact]
        public unsafe void DisplayInterceptsAreTheClosestPasses()
        {
            (CelestialBody earth, _) = Bodies.MakeEarthMoon();

            Baluev.MoidInfo* info = stackalloc Baluev.MoidInfo[6];

            int checkedPairs = 0, moreThanTwo = 0;

            foreach ((Orbit p, Orbit s) in RandomOrbitPairs(earth, 7, 500))
            {
                p.StartUT = 0;

                double Separation(double ta)
                {
                    double ut = p.StartUT + p.GetDTforTrueAnomaly(ta, 0.0);
                    return (p.getRelativePositionAtUT(ut) - s.getRelativePositionAtUT(ut)).magnitude;
                }

                int num = HarmonyPatches.FindMoidMinima(p, s, info);

                // these fall back to the stock solver
                if (num < 1)
                    continue;

                num = FindClosestPointsAdapter.SameSideMinima(p, s, info, num);

                Result r = Call(InterceptDisplay.SelectDisplayIntercepts, p, s);

                Assert.Equal(Math.Min(num, 2), r.Count);

                if (num < 1)
                    continue;

                var separations = new List<double>();
                for (int i = 0; i < num; i++)
                    separations.Add(Separation(info[i].ta1));
                separations.Sort();

                Assert.Equal(separations[0], Separation(r.FFp), 6);
                if (r.Count > 1)
                    Assert.Equal(separations[1], Separation(r.SFp), 6);

                if (num > 2)
                    moreThanTwo++;
                checkedPairs++;
            }

            _output.WriteLine($"checked {checkedPairs} orbit pairs, {moreThanTwo} with more than two minima");
            Assert.True(checkedPairs > 400);
        }

        // PatchProcessor.GetOriginalInstructions() needs Harmony's shared state, which can't be initialized in the test host, but
        // ReadMethodBody() gives the same opcodes and operands (without labels).
        private static List<CodeInstruction> ReadInstructions(MethodBase method) =>
            PatchProcessor.ReadMethodBody(method).Select(kv => new CodeInstruction(kv.Key, kv.Value)).ToList();

        [Theory]
        [InlineData(typeof(OrbitUtil), nameof(OrbitUtil.CalculateIntersections))]
        [InlineData(typeof(OrbitTargeter), "UpdateNodesAndVectors")]
        public void TranspilerRedirectsTheStockCall(Type type, string method)
        {
            MethodBase original = AccessTools.Method(type, method);
            FieldInfo  field    = AccessTools.Field(typeof(Orbit), nameof(Orbit.FindClosestPoints));
            MethodInfo invoke   = AccessTools.Method(typeof(Orbit.FindClosestPointsDelegate), "Invoke");
            MethodInfo select   = AccessTools.Method(typeof(InterceptDisplay), nameof(InterceptDisplay.SelectDisplayIntercepts));

            List<CodeInstruction> before = ReadInstructions(original);
            Assert.Single(before, c => c.LoadsField(field));

            List<CodeInstruction> after = InterceptDisplay.RedirectFindClosestPoints(ReadInstructions(original), method).ToList();

            Assert.Equal(before.Count, after.Count);
            Assert.DoesNotContain(after, c => c.LoadsField(field));
            Assert.DoesNotContain(after, c => c.Calls(invoke));
            Assert.Single(after, c => c.Calls(select));
        }

        [Fact]
        public void TranspilerLeavesUnexpectedMethodsAlone()
        {
            MethodBase original = AccessTools.Method(typeof(Orbit), nameof(Orbit._SolveSOI));

            List<CodeInstruction> before = ReadInstructions(original);
            List<CodeInstruction> after  = InterceptDisplay.RedirectFindClosestPoints(before, "Orbit._SolveSOI").ToList();

            Assert.Equal(before.Select(c => (c.opcode, c.operand)), after.Select(c => (c.opcode, c.operand)));
        }
    }
}
