using System;

namespace Csma
{
    public static class Maths
    {
        public static readonly double EPS = Math.Pow(2.0, -52.0);
        public static readonly double TINY = Math.Pow(2.0, -966.0);

        /// <summary>
        /// sqrt (a^2+b^2) without under/overflow.
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns></returns>
        public static double Hypot(double a, double b)
        {
            double r = 0;
            var absA = Math.Abs(a);
            var absB = Math.Abs(b);
            if (absA > absB)
            {
                r = b / a;
                r = absA * Math.Sqrt(1 + r * r);
            }
            else if (0 != b)
            {
                r = a / b;
                r = absB * Math.Sqrt(1 + r * r);
            }
            return r;
        }
    }
}
