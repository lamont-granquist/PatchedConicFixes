using System;
using System.Reflection;
using HarmonyLib;
using Smooth.Pools;
using Xunit;
using Xunit.Abstractions;

namespace PatchedConicFixes.Tests
{
    [Collection("KSP")]
    public class TargetingSubdivisionTests : IDisposable
    {
        private readonly Harmony _harmony = new Harmony("PatchedConicFixes.Tests.Subdivision." + Guid.NewGuid());
        private readonly ITestOutputHelper _output;
        [ThreadStatic] private static int _evaluations;
        [ThreadStatic] private static int _depth;
        [ThreadStatic] private static int _maxDepth;

        public TargetingSubdivisionTests(ITestOutputHelper output)
        {
            _output = output;
            Logger.Register(message => _output.WriteLine(message.ToString()));
            Type patch = typeof(HarmonyPatches).Assembly.GetType("PatchedConicFixes.TargetingIntervalSubdivision", true)!;
            _harmony.CreateClassProcessor(patch).Patch();
            _evaluations = _depth = _maxDepth = 0;
            _harmony.Patch(AccessTools.Method(typeof(Targeting.Sample), nameof(Targeting.Sample.Borrow)),
                prefix: new HarmonyMethod(typeof(TargetingSubdivisionTests), nameof(EvaluationBudget)));
            foreach (MethodBase method in new[] {
                AccessTools.Method(typeof(Targeting), "add_crossing_subdivisions"),
                AccessTools.Method(typeof(Targeting.Interval), nameof(Targeting.Interval.FindRoot)) })
                _harmony.Patch(method,
                    prefix: new HarmonyMethod(typeof(TargetingSubdivisionTests), nameof(RecursionBudget)),
                    finalizer: new HarmonyMethod(typeof(TargetingSubdivisionTests), nameof(LeaveRecursion)));
        }

        public void Dispose()
        {
            _harmony.UnpatchAll(_harmony.Id);
            Logger.Register(_ => { });
        }

        [Theory]
        [InlineData(3.0, -3.0)]
        [InlineData(3.1, -3.1)]
        [InlineData(-1.0, 1.0)]
        [InlineData(3.141591653589793, -3.141591653589793)]
        [InlineData(2.0, 3.0)]
        // 0.1.4 hang: repeated unchanged interval in OrbitTargeter on 10 October.
        [InlineData(1.2999800827696877, 0.099998467905360572)]
        public void SubdivisionPreservesAndShrinksAngularInterval(double start, double end)
        {
            (CelestialBody earth, _) = Bodies.MakeEarthMoon();
            var reference = new Orbit(0, 0.1, 7000000, 0, 0, 0, 0, earth);
            var target = new Orbit(3, 0.2, 7200000, 0, 20, 0, 0, earth);
            Targeting.Conic first = Targeting.Conic.Borrow(reference);
            Targeting.Conic second = Targeting.Conic.Borrow(target);
            Targeting.Interval? interval = null;
            System.Collections.Generic.List<Targeting.Interval>? children = null;
            try
            {
                Targeting.Sample left = Targeting.Sample.Borrow(first, second, start);
                Targeting.Sample right = Targeting.Sample.Borrow(first, second, end);
                Assert.True(left.valid && right.valid);
                interval = Targeting.Interval.Borrow(left, right, 0);
                children = interval.Subdivide();
                double parentWidth = ForwardWidth(start, end);
                double totalWidth = 0;
                foreach (Targeting.Interval child in children)
                {
                    double width = ForwardWidth(child.s1.v, child.s2.v);
                    Assert.InRange(width, 0, parentWidth);
                    Assert.True(width < parentWidth, "Every child must make strict progress.");
                    totalWidth += width;
                }
                Assert.Equal(parentWidth, totalWidth, 10);
            }
            finally
            {
                if (children != null)
                {
                    foreach (Targeting.Interval child in children)
                        child.Release();
                    ListPool<Targeting.Interval>.Instance.Release(children);
                }
                interval?.Release();
                Targeting.Sample.ReleaseAllBorrowed();
                first.Release();
                second.Release();
            }
        }

        private static double ForwardWidth(double start, double end)
        {
            double width = end - start;
            return width < 0 ? width + 2 * Math.PI : width;
        }

        [Fact]
        public void StartupPatchesKeepBaluevEncounterAndStockTargetingWorking()
        {
            // Exercise the same PatchAll used at startup, including the 0.1.4
            // encounter prefix and the additional stock targeting transpiler.
            _harmony.Unpatch(AccessTools.Method(typeof(Targeting.Interval), nameof(Targeting.Interval.Subdivide)),
                HarmonyPatchType.Transpiler, _harmony.Id);
            _harmony.PatchAll(typeof(HarmonyPatches).Assembly);
            _harmony.Patch(AccessTools.Method(typeof(HarmonyPatches), nameof(HarmonyPatches.CheckEncounterOldMOID)),
                prefix: new HarmonyMethod(typeof(TargetingSubdivisionTests), nameof(RejectOldEncounterSolver)));
            (CelestialBody earth, CelestialBody moon) = Bodies.MakeEarthMoon();
            var vessel = new Orbit(0.122026203095501, 0.971510969531635, 234243783.190063,
                358.94830419701, 179.012674829019, 0.000723530899600804, 432320.676460583, earth)
                { StartUT = 432320.676460583, EndUT = 1560590.81583617 };
            Assert.True(PatchedConics._CheckEncounter(vessel, new Orbit(), 432320.676460583,
                moon.orbitDriver, moon, new PatchedConics.SolverParameters(), false));
            Assert.Equal(1163413.6494660573, vessel.EndUT, 3);
            CoplanarEllipseAndCircleReturnTheAnalyticCrossings();
        }

        private static void RejectOldEncounterSolver() =>
            throw new InvalidOperationException("The 0.1.4 encounter prefix must retain the Baluev solver.");

        [Fact]
        public void CapturedIntervalViolatesSubdivisionProgressWithoutCorrection()
        {
            MethodBase subdivision = AccessTools.Method(typeof(Targeting.Interval), nameof(Targeting.Interval.Subdivide));
            Type patch = typeof(HarmonyPatches).Assembly.GetType("PatchedConicFixes.TargetingIntervalSubdivision", true)!;
            _harmony.Unpatch(subdivision, HarmonyPatchType.Transpiler, _harmony.Id);
            try
            {
                Assert.ThrowsAny<Xunit.Sdk.XunitException>(() =>
                    SubdivisionPreservesAndShrinksAngularInterval(1.2999800827696877, 0.099998467905360572));
            }
            finally { _harmony.CreateClassProcessor(patch).Patch(); }
            SubdivisionPreservesAndShrinksAngularInterval(1.2999800827696877, 0.099998467905360572);
        }

        // Test-only sentries keep a failing stock numerical search from exhausting
        // the test host's stack or looping forever. They do not ship in the mod.
        private static void EvaluationBudget()
        {
            if (++_evaluations > 20000)
                throw new InvalidOperationException("Targeting replay exceeded its sample evaluation budget.");
        }

        private static void RecursionBudget(MethodBase __originalMethod, out bool __state)
        {
            __state = false;
            if (_depth >= 128)
                throw new InvalidOperationException("Targeting replay exceeded its recursion budget in " + __originalMethod.Name + ".");
            __state = true;
            _maxDepth = Math.Max(_maxDepth, ++_depth);
        }

        private static Exception? LeaveRecursion(Exception? __exception, bool __state)
        {
            if (__state) --_depth;
            return __exception;
        }

        [Theory]
        [InlineData(0.0)]
        [InlineData(0.25)]
        [InlineData(0.5)]
        [InlineData(0.75)]
        [InlineData(1.0)]
        public void TouristACaptureNodeIntersectionsTerminate(double burnFraction)
        {
            // Quicksave #70: Tourist A and its selected Mun Alpha Station target.
            // An impulse at the saved node approximates progress through the burn;
            // it is not a replay of MechJeb's finite-burn steering or the exact hang.
            (_, CelestialBody mun) = Bodies.MakeParentChild(3.5316e12, 84159286, 6.5138398e10, 2429559,
                0, 0, 12000000, 0, 0, 0, 0);
            var before = new Orbit(0.68689060785049239, 1.2145753443843521, -1044094.1992931139,
                202.31215067973315, 82.376373426246388, -0.59012858437478566, 27639602.880920976, mun);
            var station = new Orbit(0.86493589771037149, 0.0015761190227650619, 225002.50698140604,
                62.124382992417338, 336.70175014507663, 0.53380640740365726, 27041359.646336451, mun);
            const double nodeTime = 27642069.704133388;
            Vector3d position = before.getRelativePositionAtUT(nodeTime);
            Vector3d velocity = before.getOrbitalVelocityAtUT(nodeTime);
            var during = new Orbit();
            during.UpdateFromStateVectors(position, velocity - velocity.normalized * (266.81309009934671 * burnFraction), mun, nodeTime);
            during.StartUT = nodeTime;
            try
            {
                var result = OrbitUtil.CalculateIntersections(during, station);
                _output.WriteLine($"fraction={burnFraction}; eccentricity={during.eccentricity:R}; samples={_evaluations}; maxDepth={_maxDepth}; intersections={result.numberOfIntersections}");
                Assert.InRange(result.numberOfIntersections, 0, 2);
                if (result.numberOfIntersections > 0)
                {
                    Assert.True(IsFinite(result.intersect1UT));
                    Assert.True(IsFinite(result.intersect1Distance));
                }
                if (result.numberOfIntersections > 1)
                {
                    Assert.True(IsFinite(result.intersect2UT));
                    Assert.True(IsFinite(result.intersect2Distance));
                }
            }
            finally { Targeting.Sample.ReleaseAllBorrowed(); }
        }

        private static bool IsFinite(double value) => !double.IsNaN(value) && !double.IsInfinity(value);

        [Fact]
        public void CoplanarEllipseAndCircleReturnTheAnalyticCrossings()
        {
            (CelestialBody earth, _) = Bodies.MakeEarthMoon();
            const double radius = 7000000;
            const double eccentricity = 0.2;
            var ellipse = new Orbit(0, eccentricity, radius, 0, 0, 0, 0, earth);
            var circle = new Orbit(0, 0, radius, 0, 0, 0, 0, earth);
            double distance1 = 0, distance2 = 0, anomaly1 = 0, target1 = 0, anomaly2 = 0, target2 = 0;
            int iterations = 0;
            int count = Orbit._FindClosestPoints(ellipse, circle, ref distance1, ref distance2,
                ref anomaly1, ref target1, ref anomaly2, ref target2, 0.0001, 20, ref iterations);
            Assert.Equal(2, count);
            // r=a(1-e^2)/(1+e*cos(v))=a implies cos(v)=-e.
            double expected = Math.Acos(-eccentricity);
            Assert.Equal(expected, Math.Abs(anomaly1), 7);
            Assert.Equal(expected, Math.Abs(anomaly2), 7);
            Assert.True(anomaly1 * anomaly2 < 0);
            Assert.InRange(distance1, 0, 0.1);
            Assert.InRange(distance2, 0, 0.1);
        }

        [Fact]
        public void TouristBApproachOrbitAndMunStationTerminate()
        {
            (_, CelestialBody mun) = Bodies.MakeParentChild(3.5316e12, 84159286, 6.5138398e10, 2429559,
                0, 0, 12000000, 0, 0, 0, 0);
            var vessel = new Orbit(0.74757557732379953, 1.193939495763211, -1157603.6085745883,
                5.8464171446480577, 287.60329768724006, -1.2209111188268718, 27639782.440916963, mun)
                { StartUT = 27639782.440916963 };
            var station = new Orbit(0.86493589771037149, 0.0015761190227650619, 225002.50698140604,
                62.124382992417338, 336.70175014507663, 0.53380640740365726, 27041359.646336451, mun);
            var result = OrbitUtil.CalculateIntersections(vessel, station);
            _output.WriteLine($"Tourist B: samples={_evaluations}; maxDepth={_maxDepth}; intersections={result.numberOfIntersections}");
            Assert.InRange(result.numberOfIntersections, 0, 2);
            Assert.True(_evaluations > 0);
            if (result.numberOfIntersections > 0)
                Assert.True(IsFinite(result.intersect1Distance) && IsFinite(result.intersect1UT));
        }

        [Fact]
        public void NearCircularMunStationSweep()
        {
            (_, CelestialBody mun) = Bodies.MakeParentChild(3.5316e12, 84159286, 6.5138398e10, 2429559,
                0, 0, 12000000, 0, 0, 0, 0);
            var station = new Orbit(0.86493589771037149, 0.0015761190227650619, 225002.50698140604,
                62.124382992417338, 336.70175014507663, 0.53380640740365726, 27041359.646336451, mun);
            Type patch = typeof(HarmonyPatches).Assembly.GetType("PatchedConicFixes.TargetingIntervalSubdivision", true)!;
            int stockFailures = 0;
            for (int pass = 0; pass < 2; pass++)
            {
                if (pass == 0)
                    _harmony.Unpatch(AccessTools.Method(typeof(Targeting.Interval), nameof(Targeting.Interval.Subdivide)),
                        HarmonyPatchType.Transpiler, _harmony.Id);
                else
                    _harmony.CreateClassProcessor(patch).Patch();

                var random = new Random(1970);
                for (int index = 0; index < 200; index++)
                {
                    var vessel = new Orbit(0.5 + random.NextDouble(), 0.02 * random.NextDouble(),
                        223500 + 3000 * random.NextDouble(), 360 * random.NextDouble(),
                        360 * random.NextDouble(), random.NextDouble(), 27642098.143177018, mun)
                        { StartUT = 27642098.143177018 };
                    _evaluations = _depth = _maxDepth = 0;
                    try
                    {
                        var result = OrbitUtil.CalculateIntersections(vessel, station);
                        Assert.InRange(result.numberOfIntersections, 0, 2);
                        if (result.numberOfIntersections > 0)
                            Assert.True(IsFinite(result.intersect1Distance) && IsFinite(result.intersect1UT));
                        if (result.numberOfIntersections > 1)
                            Assert.True(IsFinite(result.intersect2Distance) && IsFinite(result.intersect2UT));
                    }
                    catch (InvalidOperationException ex) when (pass == 0 && ex.Message.StartsWith("Targeting replay exceeded"))
                    {
                        stockFailures++;
                        _output.WriteLine($"stock failed case={index}; e={vessel.eccentricity:R}; a={vessel.semiMajorAxis:R}; LAN={vessel.LAN:R}; argPe={vessel.argumentOfPeriapsis:R}; inc={vessel.inclination:R}; reason={ex.Message}");
                    }
                    finally { Targeting.Sample.ReleaseAllBorrowed(); }
                }
            }
            _output.WriteLine($"Stock budget failures={stockFailures}/200; corrected completed=200/200.");
        }
    }
}
