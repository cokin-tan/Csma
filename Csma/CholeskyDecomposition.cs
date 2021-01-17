/** Cholesky Decomposition.
 * <P>
 * For a symmetric, positive definite matrix A, the Cholesky decomposition
 * is an lower triangular matrix L so that A = L*L'.
 * <P>
 * If the matrix is not symmetric or positive definite, the constructor
 * returns a partial decomposition and sets an internal flag that may
 * be queried by the isSPD() method.
 */

using System;

namespace Csma
{
    public struct CholeskyDecomposition
    {
        #region Variables
        /// <summary>
        /// Array for internal storage of decomposition.
        /// </summary>
        private double[][] _cholesky;

        /// <summary>
        /// Row and column dimension (square matrix).
        /// </summary>
        private int _n;

        /// <summary>
        /// Symmetric and positive definite flag.
        /// </summary>
        public bool isspd { get; private set; }
        #endregion

        #region Constructors
        public CholeskyDecomposition(Matrix matrix)
        {
            // Initialize
            var array = matrix.GetArray();
            this._n = matrix.row;
            this._cholesky = new double[_n][];
            for (int i = 0; i < _n; ++i)
            {
                this._cholesky[i] = new double[_n];
            }
            isspd = (matrix.column == _n);

            // Main loop
            for (int j = 0; j < _n; ++j)
            {
                var arowj = _cholesky[j];
                double d = 0.0;
                for (int k = 0; k < j; ++k)
                {
                    var arowk = _cholesky[k];
                    double s = 0.0;
                    for (int i = 0; i < k; ++i)
                    {
                        s += arowk[i] * arowj[i];
                    }
                    arowj[k] = s = (array[j][k] - s) / _cholesky[k][k];
                    d += s * s;
                    isspd = isspd & (array[k][j] == array[j][k]);
                }
                d = array[j][j] - d;
                isspd = isspd & (d > 0);
                _cholesky[j][j] = Math.Sqrt(Math.Max(d, 0));
                for (int k = j + 1; k < _n; ++k)
                {
                    _cholesky[j][k] = 0;
                }
            }
        }
        #endregion

        #region Public Methods
        /// <summary>
        /// Return triangular factor.
        /// </summary>
        /// <returns></returns>
        public Matrix GetL()
        {
            return new Matrix(_cholesky, _n, _n);
        }

        /// <summary>
        /// Solve A*X = B
        /// </summary>
        /// <param name="b">A Matrix with as many rows as A and any number of columns.</param>
        /// <returns>X so that L*L'*X = B</returns>
        /// <exception>ArgumentException : Matrix row dimensions must agree.</exception>
        /// <exception>Exception : Matrix is nor symmetric positive definite.</exception>
        public Matrix Solve(Matrix b)
        {
            if (b.row != _n)
            {
                throw new ArgumentException("Matrix row dimensions must agree.");
            }
            if (!isspd)
            {
                throw new Exception("Matrix is nor symmetric positive definite.");
            }

            // Copy right hand side.
            var x = b.GetArrayCopy();
            int nx = b.column;
            // Solve L*Y = B
            for (int k = 0; k < _n; ++k)
            {
                for (int j = 0; j < nx; ++j)
                {
                    for (int i = 0; i < k; ++i)
                    {
                        x[k][j] -= x[i][j] * _cholesky[k][i];
                    }
                    x[k][j] /= _cholesky[k][k];
                }
            }

            // Solve L'*X = Y
            for (int k = _n - 1; k >= 0; --k)
            {
                for (int j = 0; j < nx; ++j)
                {
                    for (int i = k + 1; i < _n; ++i)
                    {
                        x[k][j] -= x[i][j] * _cholesky[i][k];
                    }
                    x[k][j] /= _cholesky[k][k];
                }
            }

            return new Matrix(x, _n, nx);
        }
        #endregion
    }
}
