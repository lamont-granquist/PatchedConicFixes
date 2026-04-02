using Unity.Burst;
using Unity.Mathematics;

namespace PatchedConicFixes
{
    [BurstCompile]
    public struct COrbitData
    {
        public double a,  e,  i, w, Om;
        public double P0, P1, P2; // P vector
        public double Q0, Q1, Q2; // Q vector

        public COrbitData(double a, double e, double i, double w, double Om)
        {
            this.a = a; this.e = e; this.i = i; this.w = w; this.Om = Om;

            double cw = math.cos(w),  sw = math.sin(w);
            double ci = math.cos(i),  si = math.sin(i);
            double cO = math.cos(Om), sO = math.sin(Om);

            P0 = cw * cO - ci * sw * sO;
            P1 = cw * sO + ci * sw * cO;
            P2 = si * sw;
            Q0 = -sw * cO - ci * cw * sO;
            Q1 = ci * cw * cO - sw * sO;
            Q2 = si * cw;
        }
    }
}
