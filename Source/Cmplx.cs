using System;
using System.Runtime.CompilerServices;
using Unity.Burst;
using Unity.Mathematics;

namespace PatchedConicFixes
{
    [BurstCompile]
    public struct Cmplx
    {
        public double Re;
        public double Im;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Cmplx(double re, double im = 0.0)
        {
            Re = re;
            Im = im;
        }

        // |z|² — deliberately not called Magnitude
        public double Norm
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => Re * Re + Im * Im;
        }

        public double Abs
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => math.sqrt(Re * Re + Im * Im);
        }

        public double Phase
        {
            [MethodImpl(MethodImplOptions.AggressiveInlining)]
            get => math.atan2(Im, Re);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Cmplx Conjugate() => new Cmplx(Re, -Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public Cmplx Inverse()
        {
            double n = Norm;
            return new Cmplx(Re / n, -Im / n);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx Polar(double r, double theta)
            => new Cmplx(r * math.cos(theta), r * math.sin(theta));

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator +(Cmplx a, Cmplx b) => new Cmplx(a.Re + b.Re, a.Im + b.Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator -(Cmplx a, Cmplx b) => new Cmplx(a.Re - b.Re, a.Im - b.Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator -(Cmplx a) => new Cmplx(-a.Re, -a.Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator *(Cmplx a, Cmplx b)
            => new Cmplx(a.Re * b.Re - a.Im * b.Im,
                a.Re * b.Im + a.Im * b.Re);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator /(Cmplx a, Cmplx b)
        {
            double d = b.Re * b.Re + b.Im * b.Im;
            return new Cmplx((a.Re * b.Re + a.Im * b.Im) / d,
                (a.Im * b.Re - a.Re * b.Im) / d);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator +(Cmplx a, double b) => new Cmplx(a.Re + b, a.Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator +(double a, Cmplx b) => new Cmplx(a + b.Re, b.Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator -(Cmplx a, double b) => new Cmplx(a.Re - b, a.Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator -(double a, Cmplx b) => new Cmplx(a - b.Re, -b.Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator *(Cmplx a, double b) => new Cmplx(a.Re * b, a.Im * b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator *(double a, Cmplx b) => new Cmplx(a * b.Re, a * b.Im);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator /(Cmplx a, double b) => new Cmplx(a.Re / b, a.Im / b);

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx operator /(double a, Cmplx b)
        {
            double d = b.Re * b.Re + b.Im * b.Im;
            return new Cmplx(a * b.Re / d, -a * b.Im / d);
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(double a, Cmplx b) => b.Re == a && b.Im == 0.0;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(double a, Cmplx b) => b.Re != a || b.Im != 0.0;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Cmplx a, double b) => a.Re == b && a.Im == 0.0;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Cmplx a, double b) => a.Re != b || a.Im != 0.0;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator ==(Cmplx a, Cmplx b) => a.Re == b.Re && a.Im == b.Im;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static bool operator !=(Cmplx a, Cmplx b) => a.Re != b.Re || a.Im != b.Im;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static implicit operator Cmplx(double x) => new Cmplx(x);

        public override bool Equals(object? obj) => obj is Cmplx other && this == other;

        public override int GetHashCode()
        {
            unchecked { return (Re.GetHashCode() * 397) ^ Im.GetHashCode(); }
        }

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx Sqr(Cmplx z) => z * z;

        [MethodImpl(MethodImplOptions.AggressiveInlining)]
        public static Cmplx Sqrt(Cmplx z)
        {
            double m = z.Abs;
            if (m == 0.0) return default;
            double t = math.sqrt((m + math.abs(z.Re)) * 0.5);
            if (z.Re >= 0.0)
                return new Cmplx(t, z.Im / (2.0 * t));
            return new Cmplx(math.abs(z.Im) / (2.0 * t),
                z.Im >= 0.0 ? t : -t);
        }

        public override string ToString()
        {
            if (Im == 0.0) return $"{Re:G17}";
            if (Re == 0.0) return $"{Im:G17}i";
            return $"{Re:G17} {(Im >= 0 ? "+" : "-")} {Math.Abs(Im):G17}i";
        }
    }
}
