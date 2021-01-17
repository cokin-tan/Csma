/*
 *LU Decomposition.
 * <P>
 *     For an m-by-n matrix A with m >= n, the LU decomposition is an m-by-n
 *     unit lower triangular matrix L, an n-by-n upper triangular matrix U,
 *     and a permutation vector piv of length m so that A(piv,:) = L*U.
 *     If m < n, then L is m-by-m and U is m-by-n.
 * <P>
 *     The LU decomposition with pivoting always exists, even if the matrix is
 *     singular, so the constructor will never fail.  The primary use of the
 *     LU decomposition is in the solution of square systems of simultaneous
 *     linear equations.  This will fail if isNonsingular() returns false.
 */

using System;

namespace Csma
{
    public struct LUDecomposition
    {
        #region Variables
        /// <summary>
        /// Array for internal storage of decomposition
        /// </summary>
        private double[][] _lu;

        /// <summary>
        /// Row and column dimensions, and pivot sign.
        /// </summary>
        /// column dimension
        /// row dimension
        /// pivot sign
        private int _m, _n, _pivsign;

        /// <summary>
        /// Internal storage of pivot vector.
        /// </summary>
        private int[] _piv;
        #endregion

        #region Constructors
        /// <summary>
        /// LU Decomposition
        /// </summary>
        /// <param name="a">Rectangular matrix</param>
        public LUDecomposition(Matrix a)
        {
            // Use a "left-looking", dot-product, Crout/Doolittle algorithm.
            _lu = a.GetArrayCopy();
            _m = a.row;
            _n = a.column;
            _piv = new int[_m];
            for (int i = 0; i < _m; ++i)
            {
                _piv[i] = i;
            }
            _pivsign = 1;
            double[] lurowi;
            double[] lucolj = new double[_m];

            // Outer loop
            for (int j = 0; j < _n; ++j)
            {
                // Make a copy of the j-th column to localize references.
                for (int i = 0; i < _m; ++i)
                {
                    lucolj[i] = _lu[i][j];
                }

                // Apply previous transformations.
                for (int i = 0; i < _m; ++i)
                {
                    lurowi = _lu[i];

                    // Most of the time is spent in the following dot product.

                    int kmax = Math.Min(i, j);
                    double s = 0;
                    for (int k = 0; k < kmax; ++k)
                    {
                        s += lurowi[k] * lucolj[k];
                    }
                    lurowi[j] = lucolj[i] -= s;
                }

                // Find pivot and exchange if necessary
                int p = j;
                for (int i = j + 1; i < _m; i++)
                {
                    if (Math.Abs(lucolj[i]) > Math.Abs(lucolj[p]))
                    {
                        p = i;
                    }
                }
                if (p != j)
                {
                    for (int kk = 0; kk < _n; ++kk)
                    {
                        double t = _lu[p][kk]; 
                        _lu[p][kk] = _lu[j][kk]; 
                        _lu[j][kk] = t;
                    }
                    int k = _piv[p]; 
                    _piv[p] = _piv[j]; 
                    _piv[j] = k;
                    _pivsign = -_pivsign;
                }

                // Compute multipliers.
                if (j < _m & _lu[j][j] != 0)
                {
                    for (int i = j + 1; i < _m; ++i)
                    {
                        _lu[i][j] /= _lu[j][j];
                    }
                }
            }
        }
        #endregion

        #region Public Methods
        /// <summary>
        /// Is the matrix nonsingular?
        /// </summary>
        /// <returns>true if U, and hence A, is nonsingular.</returns>
        public bool IsNonsingular()
        {
            for (int i = 0; i < _n; ++i)
            {
                if (_lu[i][i] == 0)
                {
                    return false;
                }
            }

            return true;
        }

        /// <summary>
        /// Return lower triangular factor
        /// </summary>
        /// <returns>L</returns>
        public Matrix GetL()
        {
            Matrix x = new Matrix(_m, _n);
            var L = x.GetArray();
            for (int i = 0; i < _m; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    if (i > j)
                    {
                        L[i][j] = _lu[i][j];
                    }
                    else if (i == j)
                    {
                        L[i][j] = 1.0;
                    }
                    else
                    {
                        L[i][j] = 0.0;
                    }
                }
            }
            return x;
        }

        /// <summary>
        /// Return upper triangular factor
        /// </summary>
        /// <returns>U</returns>
        public Matrix GetU()
        {
            Matrix x = new Matrix(_n, _n);
            var U = x.GetArray();
            for (int i = 0; i < _n; ++i)
            {
                for (int j = 0; j < _n; ++j)
                {
                    if (i <= j)
                    {
                        U[i][j] = _lu[i][j];
                    }
                    else
                    {
                        U[i][j] = 0.0;
                    }
                }
            }
            return x;
        }

        /// <summary>
        /// Return pivot permutation vector
        /// </summary>
        /// <returns>piv</returns>
        public int[] GetPivot()
        {
            int[] p = new int[_m];
            for (int i = 0; i < _m; ++i)
            {
                p[i] = _piv[i];
            }
            return p;
        }

        /// <summary>
        /// Return pivot permutation vector as a one-dimensional double array
        /// </summary>
        /// <returns>double piv</returns>
        public double[] GetDoublePivot()
        {
            double[] p = new double[_m];
            for (int i = 0; i < _m; ++i)
            {
                p[i] = (double)_piv[i];
            }
            return p;
        }

        /// <summary>
        /// Determinant
        /// </summary>
        /// <returns>det(A)</returns>
        /// <exception>ArgumentException : Matrix must be square.</exception>
        public double Det()
        {
            if (_m != _n)
            {
                throw new ArgumentException("Matrix must be square.");
            }
            double d = (double)_pivsign;
            for (int j = 0; j < _n; j++)
            {
                d *= _lu[j][j];
            }
            return d;
        }

        /// <summary>
        /// Solve A*X = B
        /// </summary>
        /// <param name="b">A Matrix with as many rows as A and any number of columns.</param>
        /// <returns>X so that L*U*X = B(piv,:)</returns>
        /// <exception>ArgumentException : Matrix row dimensions must agree.</exception>
        /// <exception>Exception : Matrix is singular.</exception>
        public Matrix Solve(Matrix b)
        {
            if (b.row != _m)
            {
                throw new ArgumentException("Matrix row dimensions must agree.");
            }
            if (!this.IsNonsingular())
            {
                throw new Exception("Matrix is singular.");
            }

            // Copy right hand side with pivoting
            var nx = b.column;
            var xmat = b.GetMatrix(_piv, 0, nx - 1);
            var X = xmat.GetArray();

            // Solve L*Y = B(piv,:)
            for (int k = 0; k < _n; ++k)
            {
                for (int i = k + 1; i < _n; ++i)
                {
                    for (int j = 0; j < nx; ++j)
                    {
                        X[i][j] -= X[k][j] * _lu[i][k];
                    }
                }
            }
            // Solve U*X = Y;
            for (int k = _n - 1; k >= 0; --k)
            {
                for (int j = 0; j < nx; ++j)
                {
                    X[k][j] /= _lu[k][k];
                }
                for (int i = 0; i < k; ++i)
                {
                    for (int j = 0; j < nx; ++j)
                    {
                        X[i][j] -= X[k][j] * _lu[i][k];
                    }
                }
            }
            return xmat;
        }
        #endregion
    }
}
