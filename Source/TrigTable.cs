using Unity.Burst;
using Unity.Mathematics;

namespace PatchedConicFixes
{
    [BurstCompile]
    public unsafe struct TrigTable
    {
        public const int DIM = 21;

        public fixed double Cos[DIM];
        public fixed double Sin[DIM];

        public static TrigTable Create()
        {
            var result = new TrigTable();
            for (int j = 0; j < DIM; j++)
            {
                double angle = 2.0 * math.PI_DBL * j / DIM;
                result.Cos[j] = math.cos(angle);
                result.Sin[j] = math.sin(angle);
            }

            return result;
        }
    }
}
