/*
 *Singular Value Decomposition.
 * <P>
 *     For an m-by-n matrix A with m >= n, the singular value decomposition is
 *     an m-by-n orthogonal matrix U, an n-by-n diagonal matrix S, and
 *     an n-by-n orthogonal matrix V so that A = U*S*V'.
 * <P>
 *     The singular values, sigma[k] = S[k][k], are ordered so that
 *     sigma[0] >= sigma[1] >= ... >= sigma[n-1].
 * <P>
 *     The singular value decomposition always exists, so the constructor will
 *     never fail.  The matrix condition number and the effective numerical
 *     rank can be computed from this decomposition.
 */

using System;

namespace Csma
{
    public struct SingularValueDecomposition
    {
        #region Variables
        /// <summary>
        /// internal storage of u
        /// </summary>
        private double[][] _u;

        /// <summary>
        /// internal storage of v
        /// </summary>
        private double[][] _v;

        /// <summary>
        /// Array for internal storage of singular values.
        /// </summary>
        private double[] _s;

        /// <summary>
        /// Row dimension.
        /// </summary>
        private int _m;

        /// <summary>
        /// Column dimension.
        /// </summary>
        private int _n;
        #endregion

        #region Constructors
        /// <summary>
        /// Construct the singular value decomposition
        /// </summary>
        /// <param name="x">Rectangular matrix</param>
        public SingularValueDecomposition(Matrix x)
        {
            // Derived from LINPACK code.
            // Initialize
            var a = x.GetArrayCopy();
            _m = x.row;
            _n = x.column;

            /* Apparently the failing cases are only a proper subset of (m<n), 
             * so let's not throw error.  Correct fix to come later?
             * if (m<n) 
             * {
             * throw new ArgumentException("Csma SVD only works for m >= n"); 
             * }
             */
            int nu = Math.Min(_m, _n);
            _s = new double[Math.Min(_m + 1, _n)];
            _u = new double[_m][];
            for (int i = 0; i < _m; ++i)
            {
                _u[i] = new double[nu];
            }
            _v = new double[_n][];
            for (int i = 0; i < _n; ++i)
            {
                _v[i] = new double[_n];
            }
            var e = new double[_n];
            var work = new double[_m];
            var wantu = true;
            var wantv = true;

            // Reduce A to bidiagonal form, storing the diagonal elements
            // in s and the super-diagonal elements in e.
            var nct = Math.Min(_m - 1, _n);
            var nrt = Math.Max(0, Math.Min(_n - 2, _m));
            for (int k = 0, length = Math.Max(nct, nrt); k < Math.Max(nct, nrt); ++k)
            {
                if (k < nct)
                {
                    // Compute the transformation for the k-th column and
                    // place the k-th diagonal in s[k].
                    // Compute 2-norm of k-th column without under/overflow.
                    _s[k] = 0;
                    for (int i = k; i < _m; ++i)
                    {
                        _s[k] = Maths.Hypot(_s[k], a[i][k]);
                    }
                    if (0 != _s[k])
                    {
                        if (a[k][k] < 0)
                        {
                            _s[k] = -_s[k];
                        }
                        for (int i = k; i < _m; ++i)
                        {
                            a[i][k] /= _s[k];
                        }
                        a[k][k] += 1;
                    }
                    _s[k] = -_s[k];
                }

                for (int j = k + 1; j < _n; ++j)
                {
                    if ((k < nct) && 0 != _s[k])
                    {
                        // Apply the transformation.
                        var t = 0.0;
                        for (int i = k; i < _m; ++i)
                        {
                            t += a[i][k] * a[i][j];
                        }
                        t = -t / a[k][k];
                        for (int i = k; i < _m; ++i)
                        {
                            a[i][j] += t * a[i][k];
                        }
                    }

                    // Place the k-th row of A into e for the
                    // subsequent calculation of the row transformation.
                    e[j] = a[k][j];
                }
                if (wantu && k < nct)
                {
                    // Place the transformation in U for subsequent back
                    // multiplication.
                    for (int i = k; i < _m; ++i)
                    {
                        _u[i][k] = a[i][k];
                    }
                }

                if (k < nrt)
                {
                    // Compute the k-th row transformation and place the
                    // k-th super-diagonal in e[k].
                    // Compute 2-norm without under/overflow.
                    e[k] = 0;
                    for (int i = k + 1; i < _n; ++i)
                    {
                        e[k] = Maths.Hypot(e[k], e[i]);
                    }
                    if (0 != e[k])
                    {
                        if (e[k + 1] < 0)
                        {
                            e[k] = -e[k];
                        }
                        for (int i = k + 1; i < _n; ++i)
                        {
                            e[i] /= e[k];
                        }
                        e[k + 1] += 1;
                    }
                    e[k] = -e[k];
                    if ((k + 1 < _m) && (0 != e[k]))
                    {
                        // Apply the transformation.
                        for (int i = k + 1; i < _m; ++i)
                        {
                            work[i] = 0;
                        }
                        for (int j = k + 1; j < _n; ++j)
                        {
                            for (int i = k + 1; i < _m; ++i)
                            {
                                work[i] += e[j] * a[i][j];
                            }
                        }
                        for (int j = k + 1; j < _n; ++j)
                        {
                            var t = -e[j] / e[k + 1];
                            for (int i = k + 1; i < _m; ++i)
                            {
                                a[i][j] += t * work[i];
                            }
                        }
                    }
                    if (wantv)
                    {
                        // Place the transformation in V for subsequent
                        // back multiplication.
                        for (int i = k + 1; i < _n; ++i)
                        {
                            _v[i][k] = e[i];
                        }
                    }
                }
            }

            // Set up the final bidiagonal matrix or order p.
            var p = Math.Min(_n, _m + 1);
            if (nct < _n)
            {
                _s[nct] = a[nct][nct];
            }
            if (_m < p)
            {
                _s[p - 1] = 0;
            }
            if (nrt + 1 < p)
            {
                e[nrt] = a[nrt][p - 1];
            }
            e[p - 1] = 0;

            // If required, generate U.
            if (wantu)
            {
                for (int j = nct; j < nu; ++j)
                {
                    for (int i = 0; i < _m; ++i)
                    {
                        _u[i][j] = 0;
                    }
                    _u[j][j] = 1;
                }
                for (int k = nct - 1; k >= 0; --k)
                {
                    if (0 != _s[k])
                    {
                        for (int j = k + 1; j < nu; ++j)
                        {
                            var t = 0.0;
                            for (int i = k; i < _m; ++i)
                            {
                                t += _u[i][k] * _u[i][j];
                            }
                            t = -t / _u[k][k];
                            for (int i = k; i < _m; ++i)
                            {
                                _u[i][j] += t * _u[i][k];
                            }
                        }
                        for (int i = k; i < _m; ++i)
                        {
                            _u[i][k] = -_u[i][k];
                        }
                        _u[k][k] = 1 + _u[k][k];
                        for (int i = 0; i < k - 1; ++i)
                        {
                            _u[i][k] = 0;
                        }
                    }
                    else
                    {
                        for (int i = 0; i < _m; ++i)
                        {
                            _u[i][k] = 0;
                        }
                        _u[k][k] = 1;
                    }
                }
            }

            // If required, generate V.
            if (wantv)
            {
                for (int k = _n - 1; k >= 0; --k)
                {
                    if (k < nrt && e[k] != 0)
                    {
                        for (int j = k + 1; j < nu; ++j)
                        {
                            var t = 0.0;
                            for (int i = k + 1; i < _n; ++i)
                            {
                                t += _v[i][k] * _v[i][j];
                            }
                            t = - t / _v[k + 1][k];
                            for (int i = k + 1; i < _n; ++i)
                            {
                                _v[i][j] += t * _v[i][k];
                            }
                        }
                    }
                    for (int i = 0; i < _n; ++i)
                    {
                        _v[i][k] = 0;
                    }
                    _v[k][k] = 1;
                }
            }

            // Main iteration loop for the singular values.
            var pp = p - 1;
            var iter = 0;
            var eps = Maths.EPS;
            var tiny = Maths.TINY;
            while (p > 0)
            {
                int k = 0, kase = 0;

                // Here is where a test for too many iterations would go.

                // This section of the program inspects for
                // negligible elements in the s and e arrays.  On
                // completion the variables kase and k are set as follows.

                // kase = 1     if s(p) and e[k-1] are negligible and k<p
                // kase = 2     if s(k) is negligible and k<p
                // kase = 3     if e[k-1] is negligible, k<p, and
                //              s(k), ..., s(p) are not negligible (qr step).
                // kase = 4     if e(p-1) is negligible (convergence).
                for (k = p - 2; k >= -1; --k)
                {
                    if (k == -1)
                    {
                        break;
                    }
                    if (Math.Abs(e[k]) <= (tiny + eps*(Math.Abs(_s[k]) + Math.Abs(_s[k+1]))))
                    {
                        e[k] = 0;
                        break;
                    }
                }
                if (k == p-2)
                {
                    kase = 4;
                }
                else
                {
                    var ks = 0;
                    for (ks = p - 1; ks >= k; --ks)
                    {
                        if (ks == k)
                        {
                            break;
                        }
                        var t = (ks != p ? Math.Abs(e[ks]) : 0) +
                                (ks != k + 1 ? Math.Abs(e[ks - 1]) : 0);
                        if (Math.Abs(_s[ks]) <= (tiny + eps*t))
                        {
                            _s[ks] = 0;
                            break;
                        }
                    }
                    if (ks == k)
                    {
                        kase = 3;
                    }
                    else if (ks == p - 1)
                    {
                        kase = 1;
                    }
                    else
                    {
                        kase = 2;
                        k = ks;
                    }
                }
                ++k;

                // Perform the task indicated by kase.
                switch (kase)
                {
                    // Deflate negligible s(p).
                    case 1:
                        {
                            var f = e[p - 2];
                            e[p - 2] = 0;
                            for (int j = p - 2; j >= k; --j)
                            {
                                var t = Maths.Hypot(_s[j], f);
                                var cs = _s[j] / t;
                                var sn = f / t;
                                _s[j] = t;
                                if (j != k)
                                {
                                    f = -sn * e[j - 1];
                                    e[j - 1] = cs * e[j - 1];
                                }
                                if (wantv)
                                {
                                    for (int i = 0; i < _n; ++i)
                                    {
                                        t = cs * _v[i][j] + sn * _v[i][p - 1];
                                        _v[i][p - 1] = -sn * _v[i][j] + cs * _v[i][p - 1];
                                        _v[i][j] = t;
                                    }
                                }
                            }
                        }
                        break;

                    // Split at negligible s(k).
                    case 2:
                        {
                            var f = e[k - 1];
                            e[k - 1] = 0;
                            for (int j = k; j < p; ++j)
                            {
                                var t = Maths.Hypot(_s[j], f);
                                var cs = _s[j] / t;
                                var sn = f / t;
                                _s[j] = t;
                                f = -sn * e[j];
                                e[j] = cs * e[j];
                                if (wantu)
                                {
                                    for (int i = 0; i < _m; ++i)
                                    {
                                        t = cs * _u[i][j] + sn * _u[i][k - 1];
                                        _u[i][k - 1] = -sn * _u[i][j] + cs * _u[i][k - 1];
                                        _u[i][j] = t;
                                    }
                                }
                            }
                        }
                        break;

                    // Perform one qr step.
                    case 3:
                        {
                            var scale = Math.Max(Math.Max(Math.Max(Math.Max(
                                Math.Abs(_s[p-1]), Math.Abs(_s[p-2])), Math.Abs(e[p-2])),
                                Math.Abs(_s[k])), Math.Abs(e[k]));
                            var sp = _s[p - 1] / scale;
                            var spm1 = _s[p - 2] / scale;
                            var epm1 = e[p - 2] / scale;
                            var sk = _s[k] / scale;
                            var ek = e[k] / scale;
                            var b = ((spm1 + sp) * (spm1 - sp) + epm1 * epm1) * 0.5;
                            var c = (sp * epm1) * (sp * epm1);
                            var shift = 0.0;
                            if ((0 != b) | (0 != c))
                            {
                                shift = Math.Sqrt(b * b + c);
                                if (b < 0)
                                {
                                    shift = -shift;
                                }
                                shift = c / (b + shift);
                            }
                            var f = (sk + sp) * (sk - sp) + shift;
                            var g = sk * ek;

                            // Chase zeros.
                            for (int j = k; j < p-1; ++j)
                            {
                                var t = Maths.Hypot(f, g);
                                var cs = f / t;
                                var sn = g / t;
                                if (j != k)
                                {
                                    e[j - 1] = t;
                                }
                                f = cs * _s[j] + sn * e[j];
                                e[j] = cs * e[j] - sn * _s[j];
                                g = sn * _s[j + 1];
                                _s[j + 1] = cs * _s[j + 1];
                                if (wantv)
                                {
                                    for (int i = 0; i < _n; ++i)
                                    {
                                        t = cs * _v[i][j] + sn * _v[i][j + 1];
                                        _v[i][j + 1] = -sn * _v[i][j] + cs * _v[i][j + 1];
                                        _v[i][j] = t;
                                    }
                                }
                                t = Maths.Hypot(f, g);
                                cs = f / t;
                                sn = g / t;
                                _s[j] = t;
                                f = cs * e[j] + sn * _s[j + 1];
                                _s[j + 1] = -sn * e[j] + cs * _s[j + 1];
                                g = sn * e[j + 1];
                                e[j + 1] = cs * e[j + 1];
                                if (wantu && (j < _m-1))
                                {
                                    for (int i = 0; i < _m; ++i)
                                    {
                                        t = cs * _u[i][j] + sn * _u[i][j + 1];
                                        _u[i][j + 1] = -sn * _u[i][j] + cs * _u[i][j + 1];
                                        _u[i][j] = t;
                                    }
                                }
                            }
                            e[p - 2] = f;
                            ++iter;
                        }
                        break;

                    // Convergence.
                    case 4:
                        {
                            // Make the singular values positive.
                            if (_s[k] <= 0)
                            {
                                _s[k] = _s[k] < 0 ? -_s[k] : 0;
                                if (wantv)
                                {
                                    for (int i = 0; i <= pp; ++i)
                                    {
                                        _v[i][k] = _v[i][k];
                                    }
                                }
                            }

                            // Order the singular values.
                            while (k < pp)
                            {
                                if (_s[k] >= _s[k+1])
                                {
                                    break;
                                }
                                var t = _s[k];
                                _s[k] = _s[k + 1];
                                _s[k + 1] = t;
                                if (wantv && k < _n-1)
                                {
                                    for (int i = 0; i < _n; ++i)
                                    {
                                        t = _v[i][k + 1]; _v[i][k + 1] = _v[i][k]; _v[i][k] = t;
                                    }
                                }
                                if (wantu && k < _m-1)
                                {
                                    for (int i = 0; i < _m; ++i)
                                    {
                                        t = _u[i][k + 1];
                                        _u[i][k + 1] = _u[i][k];
                                        _u[i][k] = t;
                                    }
                                }
                                ++k;
                            }
                            iter = 0;
                            --p;
                        }
                        break;
                }
            }
        }
        #endregion

        #region Public Methods
        /// <summary>
        /// Return the left singular vectors
        /// </summary>
        /// <returns>U</returns>
        public Matrix GetU()
        {
            return new Matrix(_u, _m, Math.Min(_m + 1, _n));
        }

        /// <summary>
        /// Return the right singular vectors
        /// </summary>
        /// <returns>V</returns>
        public Matrix GetV()
        {
            return new Matrix(_v, _n, _n);
        }

        /// <summary>
        /// Return the one-dimensional array of singular values
        /// </summary>
        /// <returns>diagonal of S.</returns>
        public double[] GetSingularValues()
        {
            return _s;
        }

        /// <summary>
        /// Return the diagonal matrix of singular values
        /// </summary>
        /// <returns>S</returns>
        public Matrix GetS()
        {
            Matrix x = new Matrix(_n, _n);
            var s = x.GetArray();
            for (int i = 0; i < _n; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    s[i][j] = 0;
                }
                s[i][i] = _s[i];
            }
            return x;
        }

        /// <summary>
        /// Two Norm
        /// </summary>
        /// <returns>Max(S)</returns>
        public double Norm2()
        {
            return _s[0];
        }

        /// <summary>
        /// Two norm condition number
        /// </summary>
        /// <returns>Max(S)/Min(S)</returns>
        public double Cond()
        {
            return _s[0] / _s[Math.Min(_m, _n) - 1];
        }

        /// <summary>
        /// Effective numerical matrix rank
        /// </summary>
        /// <returns>Number of non-negligible singular values.</returns>
        public int Rank()
        {
            var tol = Math.Max(_m, _n) * _s[0] * Maths.EPS;
            int r = 0;
            for (int i = 0; i < _s.Length; ++i)
            {
                if (_s[i] > tol)
                {
                    ++r;
                }
            }
            return r;
        }
        #endregion
    }
}
