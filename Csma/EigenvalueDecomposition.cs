/*
 * Eigenvalues and eigenvectors of a real matrix. 
 * <P>
 *     If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
 *     diagonal and the eigenvector matrix V is orthogonal.
 *     I.e. A = V.times(D.times(V.transpose())) and 
 *     V.times(V.transpose()) equals the identity matrix.
 * <P>
 *     If A is not symmetric, then the eigenvalue matrix D is block diagonal
 *     with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
 *     lambda + i*mu, in 2-by-2 blocks, [lambda, mu; -mu, lambda].  The
 *     columns of V represent the eigenvectors in the sense that A*V = V*D,
 *     i.e. A.times(V) equals V.times(D).  The matrix V may be badly
 *     conditioned, or even singular, so the validity of the equation
 *     A = V*D*inverse(V) depends upon V.cond().
 */

using System;

namespace Csma
{
    public struct EigenvalueDecomposition
    {
        #region Variables
        /// <summary>
        /// Row and column dimension (square matrix)
        /// </summary>
        private int _n;

        /// <summary>
        /// Sysmetric flag
        /// </summary>
        private bool _isSymmetric;

        /// <summary>
        /// Arrays for internal storage of eigenvalues.
        /// </summary>
        private double[] _d;

        /// <summary>
        /// Arrays for internal storage of eigenvalues.
        /// </summary>
        private double[] _e;

        /// <summary>
        /// Arrays for internal storage of eigenvectors.
        /// </summary>
        private double[][] _v;

        /// <summary>
        /// Array for internal storage of nonsymmetric Hessenberg form.
        /// </summary>
        private double[][] _h;

        /// <summary>
        /// Working storage for nonsymmetric algorithm
        /// </summary>
        private double[] _ort;
        #endregion

        #region Constructors
        /// <summary>
        /// Check for symmetry, then construct the eigenvalue decomposition
        /// </summary>
        /// <param name="a">Square matrix</param>
        public EigenvalueDecomposition(Matrix x)
        {
            var a = x.GetArray();
            _n = x.column;
            _d = new double[_n];
            _e = new double[_n];
            _v = new double[_n][];
            for (int i = 0; i < _n; ++i)
            {
                _v[i] = new double[_n];
            }
            _h = null;
            _ort = null;
            _cdivi = 0;
            _cdivr = 0;

            _isSymmetric = true;
            for (int j = 0; j < _n && _isSymmetric; ++j)
            {
                for (int i = 0; (i < _n) & _isSymmetric; ++i)
                {
                    _isSymmetric = (a[i][j] == a[j][i]);
                }
            }

            if (_isSymmetric)
            {
                for (int i = 0; i < _n; ++i)
                {
                    Array.Copy(a[i], _v[i], _n);
                }

                // Tridiagonalize.
                Tred2();

                // Diagonalize.
                Tql2();
            }
            else
            {
                _ort = new double[_n];
                _h = new double[_n][];
                for (int i = 0; i < _n; ++i)
                {
                    _h[i] = new double[_n];
                    Array.Copy(a[i], _h[i], _n);
                }

                // Reduce to Hessenberg form.
                Orthes();

                // Reduce Hessenberg to real Schur form.
                Hqr2();
            }
        }
        #endregion

        #region Private Methods

        /// <summary>
        /// Symmetric Householder reduction to tridiagonal form.
        /// </summary>
        private void Tred2()
        {
            //  This is derived from the Algol procedures tred2 by
            //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
            //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.
            for (int i = 0; i < _n; ++i)
            {
                _d[i] = _v[_n - 1][i];
            }

            // Householder reduction to tridiagonal form.
            for (int i = _n - 1; i > 0; --i)
            {
                // Scale to avoid under/overflow.
                double scale = 0;
                double h = 0;
                for (int k = 0; k < i; ++k)
                {
                    scale += Math.Abs(_d[k]);
                }
                if (0 == scale)
                {
                    _e[i] = _d[i - 1];
                    for (int j = 0; j < i; ++j)
                    {
                        _d[j] = _v[i - 1][j];
                        _v[i][j] = 0;
                        _v[j][i] = 0;
                    }
                }
                else
                {
                    // Generate Householder vector.
                    var scaleRate = 1 / scale;
                    for (int k = 0; k < i; ++k)
                    {
                        _d[k] *= scaleRate;
                        h += _d[k] * _d[k];
                    }

                    var f = _d[i - 1];
                    var g = Math.Sqrt(h);
                    if (f > 0)
                    {
                        g = -g;
                    }
                    _e[i] = scale * g;
                    h -= f * g;
                    _d[i - 1] = f - g;
                    for (int j = 0; j < i; ++j)
                    {
                        _e[j] = 0;
                    }

                    // Apply similarity transformation to remaining columns.
                    for (int j = 0; j < i; ++j)
                    {
                        f = _d[j];
                        _v[j][i] = f;
                        g = _e[j] + _v[j][j] * f;
                        for (int k = j+1; k < i; ++k)
                        {
                            g += _v[k][j] * _d[k];
                            _e[k] += _v[k][j] * f;
                        }
                        _e[j] = g;
                    }
                    f = 0;
                    var hRate = 1 / h;
                    for (int j = 0; j < i; ++j)
                    {
                        _e[j] *= hRate;
                        f += _e[j] * _d[j];
                    }
                    var hh = f / (h + h);
                    for (int j = 0; j < i; ++j)
                    {
                        f = _d[j];
                        g = _e[j];
                        for (int k = j; k < i; ++k)
                        {
                            _v[k][j] -= (f * _e[k] + g * _d[k]);
                        }
                        _d[j] = _v[i - 1][j];
                        _v[i][j] = 0;
                    }
                }
                _d[i] = h;
            }

            // Accumulate transformations.
            for (int i = 0; i < _n-1; ++i)
            {
                _v[_n - 1][i] = _v[i][i];
                _v[i][i] = 0;
                var h = _d[i + 1];
                if (0 != h)
                {
                    var hRate = 1 / h;
                    for (int k = 0; k <= i; ++k)
                    {
                        _d[k] = _v[k][i + 1] * hRate;
                    }
                    for (int j = 0; j <= i; ++j)
                    {
                        double g = 0;
                        for (int k = 0; k <= i; ++k)
                        {
                            g += _v[k][i + 1] * _v[k][j];
                        }
                        for (int k = 0; k <= i; ++k)
                        {
                            _v[k][j] -= g * _d[k];
                        }
                    }
                }
                for (int k = 0; k <= i; ++k)
                {
                    _v[k][i + 1] = 0;
                }
            }
            for (int j = 0; j < _n; ++j)
            {
                _d[j] = _v[_n - 1][j];
                _v[_n - 1][j] = 0;
            }

            _v[_n - 1][_n - 1] = 0;
            _e[0] = 0;
        }

        private readonly static double EPS = Math.Pow(2.0, -52);
        /// <summary>
        /// Symmetric tridiagonal QL algorithm.
        /// </summary>
        private void Tql2()
        {
            //  This is derived from the Algol procedures tql2, by
            //  Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
            //  Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.

            for (int i = 1; i < _n; ++i)
            {
                _e[i - 1] = _e[i];
            }
            _e[_n - 1] = 0;

            var f = 0.0;
            var tst1 = 0.0;
            var eps = EPS;
            for (int l = 0; l < _n; ++l)
            {
                // Find small sub diagonal element
                tst1 = Math.Max(tst1, Math.Abs(_d[l]) + Math.Abs(_e[l]));
                var m = l;
                while (m < _n)
                {
                    if (Math.Abs(_e[m]) <= eps*tst1)
                    {
                        break;
                    }
                    ++m;
                }

                // If m == l, _d[l] is an eigenvalue.
                // otherwise, iterate.
                if (m > l)
                {
                    var iter = 0;
                    do
                    {
                        ++iter; // (Could check iteration count here.)

                        // Compute implicit shift.
                        var g = _d[l];
                        var p = (_d[l + 1] - g) / (2 * _e[l]);
                        var r = Maths.Hypot(0, 1);
                        if (p < 0)
                        {
                            r = -r;
                        }

                        _d[l] = _e[l] / (p + r);
                        _d[l + 1] = _e[l] * (p + r);
                        var dl1 = _d[l + 1];
                        var h = g - _d[l];
                        for (int i = l + 2; i < _n; ++i)
                        {
                            _d[i] -= h;
                        }
                        f += h;

                        // Implicit QL transformation.
                        p = _d[m];
                        var c = 1.0;
                        var c2 = c;
                        var c3 = c;
                        var el1 = _e[l + 1];
                        var s = 0.0;
                        var s2 = 0.0;
                        for (int i = m - 1; i >= l; --i)
                        {
                            c3 = c2;
                            c2 = c;
                            s2 = s;
                            g = c * _e[i];
                            h = c * p;
                            r = Maths.Hypot(p, _e[i]);
                            _e[i + 1] = s * r;
                            s = _e[i] / r;
                            c = p / r;
                            p = c * _d[i] - s * g;
                            _d[i + 1] = h + s * (c * g + s * _d[i]);

                            // Accumulate transformation.
                            for (int k = 0; k < _n; ++k)
                            {
                                h = _v[k][i + 1];
                                _v[k][i + 1] = s * _v[k][i] + c * h;
                                _v[k][i] = c * _v[k][i] - s * h;
                            }
                        }
                        p = -s * s2 * c3 * el1 * _e[l] / dl1;
                        _e[l] = s * p;
                        _d[l] = c * p;
                        // Check for convergence.
                    } while (Math.Abs(_e[l]) > eps * tst1);
                }
                _d[l] += f;
                _e[l] = 0;
            }

            // Sort eigenvalues and corresponding vectors.
            for (int i = 0; i < _n - 1; ++i)
            {
                var k = i;
                var p = _d[i];
                for (int j = i + 1; j < _n; ++j)
                {
                    if (_d[j] < p)
                    {
                        k = j;
                        p = _d[j];
                    }
                }

                if (k != i)
                {
                    _d[k] = _d[i];
                    _d[i] = p;
                    for (int j = 0; j < _n; ++j)
                    {
                        p = _v[j][i];
                        _v[j][i] = _v[j][k];
                        _v[j][k] = p;
                    }
                }
            }
        }

        /// <summary>
        /// Nonsymmetric reduction to Hessenberg form.
        /// </summary>
        private void Orthes()
        {
            //  This is derived from the Algol procedures orthes and ortran,
            //  by Martin and Wilkinson, Handbook for Auto. Comp.,
            //  Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutines in EISPACK.

            var low = 0;
            var high = _n - 1;
            for (int m = low + 1; m < high; ++m)
            {
                // Scale column.
                var scale = 0.0;
                for (int i = m; i <= high; ++i)
                {
                    scale += Math.Abs(_h[i][m - 1]);
                }
                if (scale != 0)
                {
                    // Compute Householder transformation.
                    var h = 0.0;
                    var scaleRate = 1 / scale;
                    for (int i = high; i >= m; --i)
                    {
                        _ort[i] = _h[i][m - 1] * scaleRate;
                        h += _ort[i] * _ort[i];
                    }
                    var g = Math.Sqrt(h);
                    if (_ort[m] > 0)
                    {
                        g = -g;
                    }
                    h -= _ort[m] * g;
                    _ort[m] -= g;

                    // Apply Householder similarity transformation
                    // H = (I-u*u'/h)*H*(I-u*u')/h)
                    for (int j = m; j < _n; ++j)
                    {
                        var f = 0.0;
                        for (int i = high; i >= m; --i)
                        {
                            f += _ort[i] * _h[i][j];
                        }
                        f /= h;
                        for (int i = m; i <= high; ++i)
                        {
                            _h[i][j] -= f * _ort[i];
                        }
                    }

                    for (int i = 0; i <= high; ++i)
                    {
                        var f = 0.0;
                        for (int j = high; j >= m; --j)
                        {
                            f += _ort[j] * _h[i][j];
                        }
                        f /= h;
                        for (int j = m; j <= high; ++j)
                        {
                            _h[i][j] -= f * _ort[j];
                        }
                    }
                    _ort[m] = scale * _ort[m];
                    _h[m][m - 1] = scale * g;
                }
            }

            // Accumulate transformations (Algol's ortran).
            for (int i = 0; i < _n; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    _v[i][j] = (i == j ? 1.0 : 0.0);
                }
            }

            for (int m = high - 1; m >= low + 1; --m)
            {
                if (_h[m][m-1] != 0)
                {
                    for (int i = m + 1; i <= high; ++i)
                    {
                        _ort[i] = _h[i][m - 1];
                    }
                    for (int j = m; j <= high; ++j)
                    {
                        var g = 0.0;
                        for (int i = m; i <= high; ++i)
                        {
                            g += _ort[i] * _v[i][j];
                        }
                        // Double division avoids possible underflow
                        g = (g / _ort[m]) / _h[m][m - 1];
                        for (int i = m; i <= high; ++i)
                        {
                            _v[i][j] += g * _ort[i];
                        }
                    }
                }
            }
        }

        private double _cdivr;
        private double _cdivi;
        /// <summary>
        /// Complex scalar division
        /// </summary>
        /// <param name="xr"></param>
        /// <param name="xi"></param>
        /// <param name="yr"></param>
        /// <param name="yi"></param>
        private void Cdiv(double xr, double xi, double yr, double yi)
        {
            var r = 0.0;
            var d = 0.0;
            if (Math.Abs(yr) > Math.Abs(yi))
            {
                r = yi / yr;
                d = 1 / (yr + r * yi);
                _cdivr = (xr + r * xi) * d;
                _cdivi = (xi - r * xr) * d;
            }
            else
            {
                r = yr / yi;
                d = 1 / (yi + r * yr);
                _cdivr = (r * xr + xi) * d;
                _cdivi = (r * xi - xr) * d;
            }
        }

        /// <summary>
        /// Nonsymmetric reduction from Hessenberg to real Schur form.
        /// </summary>
        private void Hqr2()
        {
            //  This is derived from the Algol procedure hqr2,
            //  by Martin and Wilkinson, Handbook for Auto. Comp.,
            //  Vol.ii-Linear Algebra, and the corresponding
            //  Fortran subroutine in EISPACK.

            // Initialize
            var nn = _n;
            var n = nn - 1;
            var low = 0;
            var high = nn - 1;
            var eps = EPS;
            var exshift = 0.0;
            double p = 0, q = 0, r = 0, s = 0, z = 0, t = 0, w = 0, x = 0, y = 0;

            // Store roots isolated by balance and compute matrix norm
            var norm = 0.0;
            for (int i = 0; i < nn; ++i)
            {
                if (i < low | i > high)
                {
                    _d[i] = _h[i][i];
                    _e[i] = 0;
                }
                for (int j = Math.Max(i - 1, 0); j < nn; ++j)
                {
                    norm += Math.Abs(_h[i][j]);
                }
            }

            // Outer loop over eigenvalue index
            var iter = 0;
            while (n >= low)
            {
                // Look for single small sub-diagonal element
                var l = n;
                while (l > low)
                {
                    s = Math.Abs(_h[l - 1][l - 1]) + Math.Abs(_h[l][l]);
                    if (s == 0)
                    {
                        s = norm;
                    }
                    if (Math.Abs(_h[l][l-1]) < eps * s)
                    {
                        break;
                    }
                    --l;
                }

                // Check for convergence
                if (l == n) // One root found
                {
                    _h[n][n] += exshift;
                    _d[n] = _h[n][n];
                    _e[n] = 0;
                    --n;
                    iter = 0;
                }
                else if (l == n-1) // Two roots found
                {
                    w = _h[n][n - 1] * _h[n - 1][n];
                    p = (_h[n - 1][n - 1] - _h[n][n]) * 0.5;
                    q = p * p + w;
                    z = Math.Sqrt(Math.Abs(q));
                    _h[n][n] += exshift;
                    _h[n - 1][n - 1] += exshift;
                    x = _h[n][n];

                    
                    if (q >= 0) // Real pair
                    {
                        if (p >= 0)
                        {
                            z = p + z;
                        }
                        else
                        {
                            z = p - z;
                        }
                        _d[n - 1] = x + z;
                        _d[n] = _d[n - 1];
                        if (0 != z)
                        {
                            _d[n] = x - w / z;
                        }
                        _e[n - 1] = 0;
                        _e[n] = 0;
                        x = _h[n][n - 1];
                        s = Math.Abs(x) + Math.Abs(z);
                        p = x / s;
                        q = z / s;
                        r = Math.Sqrt(p * p + q * q);
                        p /= r;
                        q /= r;

                        // Row modification
                        for (int j = n - 1; j < nn; ++j)
                        {
                            z = _h[n - 1][j];
                            _h[n - 1][j] = q * z + p * _h[n][j];
                            _h[n][j] = q * _h[n][j] - p * z;
                        }

                        // Column modification
                        for (int i = 0; i <= n; ++i)
                        {
                            z = _h[i][n - 1];
                            _h[i][n - 1] = q * z + p * _h[i][n];
                            _h[i][n] = q * _h[i][n] - p * z;
                        }

                        // Accumulate transformations
                        for (int i = low; i <= high; ++i)
                        {
                            z = _v[i][n - 1];
                            _v[i][n - 1] = q * z + p * _v[i][n];
                            _v[i][n] = q * _v[i][n] - p * z;
                        }
                    }
                    else // Complex pair
                    {
                        _d[n - 1] = x + p;
                        _d[n] = x + p;
                        _e[n - 1] = z;
                        _e[n] = -z;
                    }
                    n -= 2;
                    iter = 0;
                }
                else // No convergence yet
                {
                    // From shift
                    x = _h[n][n];
                    y = 0;
                    w = 0;
                    if (l < n)
                    {
                        y = _h[n - 1][n - 1];
                        w = _h[n][n - 1] * _h[n - 1][n];
                    }


                    // Wilkinson's original ad hoc shift
                    if (10 == iter) 
                    {
                        exshift += x;
                        for (int i = low; i <= n; ++i)
                        {
                            _h[i][i] -= x;
                        }
                        s = Math.Abs(_h[n][n - 1]) + Math.Abs(_h[n - 1][n - 2]);
                        x = y = 0.75 * s;
                        w = -0.4375 * s * s;
                    }

                    // MATLAB's new ad hoc shift
                    if (30 == iter) 
                    {
                        s = (y - x) * 0.5;
                        s = s * s + w;
                        if (s > 0)
                        {
                            s = Math.Sqrt(s);
                            if (y < x)
                            {
                                s = -s;
                            }
                            s = x - w / ((y - x) * 0.5 + s);
                            for (int i = low; i <= n; ++i)
                            {
                                _h[i][i] -= s;
                            }
                            exshift += s;
                            x = y = w = 0.964;
                        }
                    }

                    ++iter; // (Could check iteration count here.)

                    // Look for two consecutive small sub-diagonal elements
                    var m = n - 2;
                    while (m >= l)
                    {
                        z = _h[m][m];
                        r = x - z;
                        s = y - z;
                        p = (r * s - w) / _h[m + 1][m] + _h[m][m + 1];
                        q = _h[m + 1][m + 1] - z - r - s;
                        r = _h[m + 2][m + 1];
                        s = Math.Abs(p) + Math.Abs(q) + Math.Abs(r);
                        p /= s;
                        q /= s;
                        r /= s;
                        if (l == m)
                        {
                            break;
                        }
                        if (Math.Abs(_h[m][m-1]) * (Math.Abs(q) + Math.Abs(r)) < 
                            eps * (Math.Abs(p) * (Math.Abs(_h[m-1][m-1]) + Math.Abs(z) +
                            Math.Abs(_h[m + 1][m + 1]))))
                        {
                            break;
                        }
                        --m;
                    }

                    for (int i = m+2; i <= n; ++i)
                    {
                        _h[i][i - 2] = 0;
                        if (i > m+2)
                        {
                            _h[i][i - 3] = 0;
                        }
                    }

                    // Double QR step involving rows l:n and columns m:n
                    for (int k = m; k <= n-1; ++k)
                    {
                        var notlast = k != n - 1;
                        if (k != m)
                        {
                            p = _h[k][k - 1];
                            q = _h[k + 1][k - 1];
                            r = notlast ? _h[k + 2][k - 1] : 0;
                            x = Math.Abs(p) + Math.Abs(q) + Math.Abs(r);
                            if (0 != x)
                            {
                                p /= x;
                                q /= x;
                                r /= x;
                            }
                        }
                        if (0 == x)
                        {
                            break;
                        }

                        s = Math.Sqrt(p * p + q * q + r * r);
                        if (p < 0)
                        {
                            s = -s;
                        }
                        if (0 != s)
                        {
                            if (k != m)
                            {
                                _h[k][k - 1] = -s * x;
                            }
                            else if (l != m)
                            {
                                _h[k][k - 1] = -_h[k][k - 1];
                            }

                            p += s;
                            x = p / s;
                            y = q / s;
                            z = r / s;
                            q /= p;
                            r /= p;

                            // Row modification
                            for (int j = k; j < nn; ++j)
                            {
                                p = _h[k][j] + q * _h[k + 1][j];
                                if (notlast)
                                {
                                    p += r * _h[k + 2][j];
                                    _h[k + 2][j] -= p * z;
                                }
                                _h[k][j] -= p * x;
                                _h[k + 1][j] -= p * y;
                            }

                            // Column modification
                            for (int i = 0, length = Math.Min(n, k + 3); i <= length; ++i)
                            {
                                p = x * _h[i][k] + y * _h[i][k + 1];
                                if (notlast)
                                {
                                    p += z * _h[i][k + 2];
                                    _h[i][k + 2] -= p * r;
                                }
                                _h[i][k] -= p;
                                _h[i][k + 1] -= p * q;
                            }

                            // Accumulate transformations
                            for (int i = low; i <= high; ++i)
                            {
                                p = x * _v[i][k] + y * _v[i][k + 1];
                                if (notlast)
                                {
                                    p += z * _v[i][k + 2];
                                    _v[i][k + 2] -= p * r;
                                }
                                _v[i][k] -= p;
                                _v[i][k + 1] -= p * q;
                            }
                        } // (0 != s)
                    }// k loop
                } // check convergence
            } // while (n >= low)

            // Backsubstitute to find vectors of upper triangular form
            if (0 == norm)
            {
                return;
            }

            for (n = nn - 1; n >= 0; --n)
            {
                p = _d[n];
                q = _e[n];

                if (0 == q) // Real vector
                {
                    var l = n;
                    _h[n][n] = 1;
                    for (int i = n - 1; i >= 0; --i)
                    {
                        w = _h[i][i] - p;
                        r = 0;
                        for (int j = l; j <= n; ++j)
                        {
                            r += _h[i][j] * _h[j][n];
                        }
                        if (_e[i] < 0)
                        {
                            z = w;
                            s = r;
                        }
                        else
                        {
                            l = i;
                            if (0 == _e[i])
                            {
                                if (0 != w)
                                {
                                    _h[i][n] = -r / w;
                                }
                                else
                                {
                                    _h[i][n] = -r / (eps * norm);
                                }
                            }
                            else // Solve real equations
                            {
                                x = _h[i][i + 1];
                                y = _h[i + 1][i];
                                q = (_d[i] - p) * (_d[i] - p) + _e[i] * _e[i];
                                t = (x * s - z * r) / q;
                                _h[i][n] = t;
                                if (Math.Abs(x) > Math.Abs(z))
                                {
                                    _h[i + 1][n] = (-r - w * t) / x;
                                }
                                else
                                {
                                    _h[i + 1][n] = (-s - y * t) / z;
                                }
                            }

                            // Overflow Control
                            t = Math.Abs(_h[i][n]);
                            if((eps * t) * t > 1)
                            {
                                for (int j = i; j <= n; ++j)
                                {
                                    _h[j][n] /= t;
                                }
                            }
                        }
                    }
                }
                else if (q < 0) // Complex vector
                {
                    var l = n - 1;

                    // Last vector component imaginary so matrix is triangular
                    if (Math.Abs(_h[n][n-1]) > Math.Abs(_h[n-1][n]))
                    {
                        _h[n - 1][n - 1] = q / _h[n][n - 1];
                        _h[n - 1][n] = (p - _h[n][n]) / _h[n][n - 1];
                    }
                    else
                    {
                        Cdiv(0, -_h[n - 1][n], _h[n - 1][n - 1] - p, q);
                        _h[n - 1][n - 1] = _cdivr;
                        _h[n - 1][n] = _cdivi;
                    }
                    _h[n][n - 1] = 0;
                    _h[n][n] = 1;
                    for (int i = n-2; i >= 0; --i)
                    {
                        double ra = 0, sa = 0, vr = 0, vi = 0;
                        for (int j = l; j <= n; ++j)
                        {
                            ra += _h[i][j] * _h[j][n - 1];
                            sa += _h[i][j] * _h[j][n];
                        }
                        w = _h[i][i] - p;

                        if (_e[i] < 0)
                        {
                            z = w;
                            r = ra;
                            s = sa;
                        }
                        else
                        {
                            l = i;
                            if (0 == _e[i])
                            {
                                Cdiv(-ra, -sa, w, q);
                                _h[i][n - 1] = _cdivr;
                                _h[i][n] = _cdivi;
                            }
                            else
                            {
                                // Solve complex equations
                                x = _h[i][i + 1];
                                y = _h[i + 1][i];
                                vr = (_d[i] - p) * (_d[i] - p) + _e[i] * _e[i] - q * q;
                                vi = (_d[i] - p) * 2 * q;
                                if (vr == 0 & vi == 0)
                                {
                                    vr = eps * norm * (Math.Abs(w) + Math.Abs(q) +
                                        Math.Abs(x) + Math.Abs(y) + Math.Abs(z));
                                }
                                Cdiv(x * r - z * ra + q * sa, x * s - z * sa - q * ra, vr, vi);
                                _h[i][n - 1] = _cdivr;
                                _h[i][n] = _cdivi;
                                if (Math.Abs(x) > (Math.Abs(z) + Math.Abs(q)))
                                {
                                    _h[i + 1][n - 1] = (-ra - w * _h[i][n - 1] + q * _h[i][n]) / x;
                                    _h[i + 1][n] = (-sa - w * _h[i][n] - q * _h[i][n - 1]) / x;
                                }
                                else
                                {
                                    Cdiv(-r - y * _h[i][n - 1], -s - y * _h[i][n], z, q);
                                    _h[i + 1][n - 1] = _cdivr;
                                    _h[i + 1][n] = _cdivi;
                                }
                            }

                            // Overflow control
                            t = Math.Max(Math.Abs(_h[i][n - 1]), Math.Abs(_h[i][n]));
                            if ((eps*t)*t > 1)
                            {
                                for (int j = i; j <= n; ++j)
                                {
                                    _h[j][n - 1] /= t;
                                    _h[j][n] /= t;
                                }
                            }
                        }
                    }
                }
            }

            // Vectors of isolated roots
            for (int i = 0; i < nn; ++i)
            {
                if (i < low | i > high)
                {
                    for (int j = i; j < nn; ++j)
                    {
                        _v[i][j] = _h[i][j];
                    }
                }
            }

            // Back transformation to get eigenvectors of original matrix
            for (int j = nn - 1; j >= low; --j)
            {
                for (int i = low; i <= high; ++i)
                {
                    z = 0;
                    for (int k = low; k <= Math.Min(j, high); ++k)
                    {
                        z += _v[i][k] * _h[k][j];
                    }
                    _v[i][j] = z;
                }
            }
        }
        #endregion

        #region Public Methods
        public Matrix GetV()
        {
            return new Matrix(_v, _n, _n);
        }

        /// <summary>
        /// Return the real parts of the eigenvalues
        /// </summary>
        /// <returns>real(diag(D))</returns>
        public double[] GetRealEigenvalues()
        {
            return _d;
        }

        /// <summary>
        /// Return the imaginary parts of the eigenvalues
        /// </summary>
        /// <returns>imag(diag(D))</returns>
        public double[] GetImagEigenvalues()
        {
            return _e;
        }

        /// <summary>
        /// Return the block diagonal eigenvalue matrix
        /// </summary>
        /// <returns>D</returns>
        public Matrix GetD()
        {
            var x = new Matrix(_n, _n);
            var d = x.GetArray();
            for (int i = 0; i < _n; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    d[i][j] = 0;
                }
                d[i][i] = _d[i];
                if (_e[i] > 0)
                {
                    d[i][i + 1] = _e[i];
                }
                else if (_e[i] < 0)
                {
                    d[i][i - 1] = _e[i];
                }
            }
            return x;
        }
        #endregion
    }
}
