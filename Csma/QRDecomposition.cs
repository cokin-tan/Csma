/*
 *QR Decomposition.
 * <P>
 *     For an m-by-n matrix A with m >= n, the QR decomposition is an m-by-n 
 *     orthogonal matrix Q and an n-by-n upper triangular matrix R so that
 *     A = Q*R.
 * <P>
 *     The QR decomposition always exists, even if the matrix does not have
 *     full rank, so the constructor will never fail.  The primary use of the
 *     QR decomposition is in the least squares solution of non-square systems
 *     of simultaneous linear equations.  This will fail if isFullRank()
 *     returns false.
 */

using System;

namespace Csma
{
    public struct QRDecomposition
    {
        /// <summary>
        /// Array for internal storage of decomposition
        /// </summary>
        private double[][] _qr;

        /// <summary>
        /// Column dimension
        /// </summary>
        private int _m;

        /// <summary>
        /// Row dimension
        /// </summary>
        private int _n;

        /// <summary>
        /// Array for internal storage of diagonal of R
        /// </summary>
        private double[] _rdiag;

        #region Constructors
        /// <summary>
        /// QR Decomposition, computed by Householder reflections.
        /// </summary>
        /// <param name="a">Rectangular matrix</param>
        public QRDecomposition(Matrix a)
        {
            this._qr = a.GetArrayCopy();
            this._m = a.row;
            this._n = a.column;
            this._rdiag = new double[_n];

            // Main loop.
            for (int k = 0; k < _n; ++k)
            {
                // Compute 2-norm of k-th column without under/overflow.
                double nrm = 0;
                for (int i = k; i < _m; ++i)
                {
                    nrm = Maths.Hypot(nrm, _qr[i][k]);
                }

                if (nrm != 0)
                {
                    // Form k-th Householder vector.
                    if (_qr[k][k] < 0)
                    {
                        nrm = -nrm;
                    }

                    for (int i = k; i < _m; ++i)
                    {
                        _qr[i][k] /= nrm;
                    }
                    _qr[k][k] += 1;

                    // Apply transformation to remaining columns.
                    for (int j = k + 1; j < _n; ++j)
                    {
                        double s = 0;
                        for (int i = k; i < _m; ++i)
                        {
                            s += _qr[i][k] * _qr[i][j];
                        }
                        s = -s / _qr[k][k];
                        for (int i = k; i < _m; ++i)
                        {
                            _qr[i][j] += s * _qr[i][k];
                        }
                    }
                }
                _rdiag[k] = -nrm;
            }
        }
        #endregion

        #region Public Methods
        /// <summary>
        /// Is the matrix full matrix?
        /// </summary>
        /// <returns>true if R, and hence A, has full rank.</returns>
        public bool IsFullRank()
        {
            for (int i = 0; i < _n; ++i)
            {
                if (0 == _rdiag[i])
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Return the Householder vectors
        /// </summary>
        /// <returns>Lower trapezoidal matrix whose columns define the reflections</returns>
        public Matrix GetH()
        {
            Matrix x = new Matrix(_m, _n);
            var h = x.GetArray();
            for (int i = 0; i < _m; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    if (i >= j)
                    {
                        h[i][j] = _qr[i][j];
                    }
                    else
                    {
                        h[i][j] = 0;
                    }
                }
            }
            return x;
        }

        /// <summary>
        /// Return the upper triangular factor
        /// </summary>
        /// <returns>R</returns>
        public Matrix GetR()
        {
            Matrix x = new Matrix(_n, _n);
            var r = x.GetArray();
            for (int i = 0; i < _n; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    if (i < j)
                    {
                        r[i][j] = _qr[i][j];
                    }
                    else if (i == j)
                    {
                        r[i][j] = _rdiag[i];
                    }
                    else
                    {
                        r[i][j] = 0;
                    }
                }
            }
            return x;
        }

        /// <summary>
        /// Generate and return the (economy-sized) orthogonal factor
        /// </summary>
        /// <returns>Q</returns>
        public Matrix GetQ()
        {
            Matrix x = new Matrix(_m, _n);
            var q = x.GetArray();
            for (int k = _n - 1; k >= 0; --k)
            {
                for (int i = 0; i < _m; ++i)
                {
                    q[i][k] = 0;
                }

                q[k][k] = 1;
                for (int j = k; j < _n; ++j)
                {
                    if (_qr[k][k] != 0)
                    {
                        double s = 0;
                        for (int i = k; i < _m; ++i)
                        {
                            s += _qr[i][k] * q[i][j];
                        }
                        s = -s / _qr[k][k];
                        for (int i = k; i < _m; ++i)
                        {
                            q[i][j] += s * _qr[i][k];
                        }
                    }
                }
            }
            return x;
        }

        /// <summary>
        /// Least squares solution of A*X = B
        /// </summary>
        /// <param name="b">A Matrix with as many rows as A and any number of columns.</param>
        /// <returns>X that minimizes the two norm of Q*R*X-B.</returns>
        /// <exception>ArgumentException : Matrix row dimensions must agree.</exception>
        /// <exception>Exception : Matrix is rank deficient.</exception>
        public Matrix Solve(Matrix b)
        {
            if (b.row != _m)
            {
                throw new ArgumentException("Matrix row dimensions must agree.");
            }
            if (!this.IsFullRank())
            {
                throw new Exception("Matrix is rank deficient.");
            }

            // Copy right hand side
            int nx = b.column;
            var x = b.GetArrayCopy();

            // Compute Y = transpose(Q)*B
            for (int k = 0; k < _n; ++k)
            {
                for (int j = 0; j < nx; ++j)
                {
                    double s = 0;
                    for (int i = k; i < _m; ++i)
                    {
                        s += _qr[i][k] * x[i][j];
                    }
                    s = -s / _qr[k][k];
                    for (int i = k; i < _m; ++i)
                    {
                        x[i][j] += s * _qr[i][k];
                    }
                }
            }

            // Solve R*X = Y
            for (int k = _n - 1; k >= 0; --k)
            {
                for (int j = 0; j < nx; ++j)
                {
                    x[k][j] /= _rdiag[k];
                }

                for (int i = 0; i < k; ++i)
                {
                    for (int j = 0; j < nx; ++j)
                    {
                        x[i][j] -= x[k][j] * _qr[i][k];
                    }
                }
            }

            return (new Matrix(x, _n, nx).GetMatrix(0, _n - 1, 0, nx - 1));
        }
        #endregion
    }
}
