/*
 * Csma = CSharp Matrix class.
 * <P>
 *     The CSharp Matrix Class provides the fundamental operations of numerical
 *     linear algebra. Various constructors create Matrices from two dimensional
 *     arrays of double precision floating point numbers. Various "gets" and 
 *     "sets" provide access to sub matrices and matrix elements. Several methods
 *     implement basic matrix arithmetic, including matrix addition and
 *     multiplication, matrix norms, and element-by-element array operations.
 *     All the operations in this version of the Matrix Class involve real matrices.
 *     Complex matrices may be handled in a future version.
 * <P>
 *     Five fundamental matrix decompositions, which consist of pairs or triples
 *     of matrices, permutation vectors, and the like, produce results in five
 *     decomposition classes. These decompositions are accessed by the Matrix
 *     class to compute solutions of simultaneous linear equations, determinants,
 *     inverses and other matrix functions. The five decompositions are:
 * <P>
 * <UL>
 *    <LI>Cholesky Decomposition of symmetric, positive definite matrices.
 *    <LI>LU Decomposition of rectangular matrices.
 *    <LI>QR Decomposition of rectangular matrices.
 *    <LI>Singular Value Decomposition of retangular matrices.
 *    <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
 * </UL>
 * <DL>
 * <DT><B>Example of use:</B></DT>
 * <P>
 * <DD>Solve a linear system A x = b and compute the residual norm,||b - A x||.
 * <P><PRE>
 *       double[][] vals = {{1,2,3},{4,5,6},{7,8,9}};
 *       Matrix A = new Matrix(vals);
 *       Matrix b = Matrix.Random(3, 1);
 *       Matrix x = A.Solve(b);
 *       Matrix r = A.Times(x) - b;
 *       double rnorm = r.NormInf();
 * </PRE></DD>
 * </DL>
 */

using System;

namespace Csma
{
    public struct Matrix
    {
        #region Variables
        private double[][] _array;

        /// <summary>
        /// Get row dimension
        /// </summary>
        /// <returns>the number of rows.</returns>
        public int row { get; private set; }

        /// <summary>
        /// Get column dimension
        /// </summary>
        /// <returns>the number of columns.</returns>
        public int column { get; private set; }
        #endregion

        #region Constructors
        /// <summary>
        /// Construct an m*n matrix of zeros
        /// </summary>
        /// <param name="m">Number of rows</param>
        /// <param name="n">Number of colums</param>
        public Matrix(int m, int n)
        {
            this.row = m;
            this.column = n;
            this._array = new double[m][];
            for (int index = 0; index < m; ++index)
            {
                this._array[index] = new double[n];
            }
        }

        /// <summary>
        /// Construct an m*n matrix of value
        /// </summary>
        /// <param name="m">Number of rows</param>
        /// <param name="n">Number of colums</param>
        /// <param name="value">Fill the matrix with this scalar value</param>
        public Matrix(int m, int n, double value)
        {
            this.row = m;
            this.column = n;
            this._array = new double[m][];
            for (int i = 0; i < m; ++i)
            {
                this._array[i] = new double[n];
                for (int j = 0; j < n; ++j)
                {
                    this._array[i][j] = value;
                }
            }
        }

        /// <summary>
        /// Construct a matrix from a 2-D array
        /// </summary>
        /// <param name="array">Two-dimensional array of double</param>
        /// <exception>ArgumentException : All rows must have the same length.</exception>
        /// <see>ConstructWithCopy</see>
        public Matrix(double[][] array)
        {
            this.row = array.Length;
            this.column = array[0].Length;
            for (int i = 0; i < this.row; ++i)
            {
                if (array[i].Length != column)
                {
                    throw new ArgumentException("All rows must have the same length");
                }
            }

            this._array = array;
        }

        /// <summary>
        /// Construct a matrix quickly without checking arguments
        /// </summary>
        /// <param name="array">Two-dimensional array of double</param>
        /// <param name="m">Number of rows</param>
        /// <param name="n">Number of colums</param>
        public Matrix(double[][] array, int m, int n)
        {
            this.row = m;
            this.column = n;
            this._array = array;
        }

        /// <summary>
        /// Construct a matrix from a one-dimensional packed array
        /// </summary>
        /// <param name="array">One-dimensional array of double, packed by columns</param>
        /// <param name="m">Number of rows</param>
        /// <exception>ArgumentException : Array length must be a multiple of m</exception>
        public Matrix(double[] array, int m)
        {
            var n = (0 != m ? array.Length / m : 0);
            if (m*n != array.Length)
            {
                throw new ArgumentException("Array length must be a multiple of m.");
            }
            this.row = m;
            this.column = n;
            this._array = new double[m][];
            for (int i = 0; i < m; ++i)
            {
                this._array[i] = new double[n];
                for (int j = 0; j < n; ++j)
                {
                    this._array[i][j] = array[i + j * m];
                }
            }
        }
        #endregion

        #region Public methods
        /// <summary>
        /// Make a deep copy of matrix
        /// </summary>
        public Matrix Copy()
        {
            var matrix = new Matrix(this.row, this.column);
            var array = matrix.GetArray();
            for (int i = 0; i < this.row; i++)
            {
                Array.Copy(_array[i], array[i], this.column);
            }

            return matrix;
        }

        /// <summary>
        /// Clone the Matrix object
        /// </summary>
        public object Clone()
        {
            return this.Copy();
        }

        /// <summary>
        /// Access the internal tow-dimensional array.
        /// </summary>
        /// <returns>Pointer to the two-dimensional array of matrix elements.</returns>
        public double[][] GetArray()
        {
            return this._array;
        }

        /// <summary>
        /// Copy the internal two-dimensional array.
        /// </summary>
        /// <returns>Two-dimensional array copy of matrix elements.</returns>
        public double[][] GetArrayCopy()
        {
            var array = new double[this.row][];
            for (int i = 0; i < this.row; ++i)
            {
                array[i] = new double[this.column];
                Array.Copy(this._array[i], array[i], this.column);
            }
            return array;
        }

        /// <summary>
        /// Make a one-dimensional column packed copy of the internal array.
        /// </summary>
        /// <returns>Matrix elements packed in a one-dimensional array by column.</returns>
        public double[] GetColumnPackedCopy()
        {
            var values = new double[row * column];
            for (int i = 0; i < this.row; ++i)
            {
                for (int j = 0; j < this.column; ++j)
                {
                    values[i + j * this.row] = _array[i][j];
                }
            }

            return values;
        }

        /// <summary>
        /// Make a one-dimensional row packed copy of the internal array.
        /// </summary>
        /// <returns>Matrix elements packed in a one-dimensional array by row.</returns>
        public double[] GetRowPackedCopy()
        {
            var values = new double[row * column];
            for (int i = 0; i < this.row; ++i)
            {
                for (int j = 0; j < this.column; ++j)
                {
                    values[i * this.column + j] = _array[i][j];
                }
            }

            return values;
        }

        /// <summary>
        /// Get a single element.
        /// </summary>
        /// <param name="i">Row index</param>
        /// <param name="j">Column index</param>
        /// <returns>_array(i,j)</returns>
        /// <exception>ArgumentOutOfRangeException</exception>
        public double Get(int i, int j)
        {
            return this._array[i][j];
        }

        /// <summary>
        /// Set a single element.
        /// </summary>
        /// <param name="i">Row index</param>
        /// <param name="j">Column index</param>
        /// <param name="value">_array(i,j) to set</param>
        /// <exception>ArgumentOutOfRangeException</exception>
        public void Set(int i, int j, double value)
        {
            _array[i][j] = value;
        }

        /// <summary>
        /// Get a sub matrix
        /// </summary>
        /// <param name="i0">Start row index</param>
        /// <param name="i1">End row index</param>
        /// <param name="j0">Start column index</param>
        /// <param name="j1">End column index</param>
        /// <returns>_array(i0:i1, j0:j1)</returns>
        /// <exception>ArgumentOutOfRangeException : sub matrix indices.</exception>
        public Matrix GetMatrix(int i0, int i1, int j0, int j1)
        {
            var matrix = new Matrix(i1 - i0 + 1, j1 - j0 + 1);
            var array = matrix.GetArray();
            try
            {
                for (int i = i0; i <= i1; ++i)
                {
                    for (int j = j0; j <= j1; ++j)
                    {
                        array[i - i0][j - j0] = _array[i][j];
                    }
                }
            }
            catch(ArgumentOutOfRangeException)
            {
                throw new ArgumentOutOfRangeException("Sub matrix indices.");
            }

            return matrix;
        }

        /// <summary>
        /// Get a sub matrix
        /// </summary>
        /// <param name="r">Array of row indices</param>
        /// <param name="c">Array of column indices</param>
        /// <returns>_array(r(:), c(:))</returns>
        /// <exception>ArgumentOutOfRangeException : sub matrix indices.</exception>
        public Matrix GetMatrix(int[] r, int[] c)
        {
            var m = r.Length;
            var n = c.Length;
            var matrix = new Matrix(m, n);
            var array = matrix.GetArray();
            try
            {
                for (int i = 0; i < m; ++i)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        array[i][j] = _array[r[i]][c[j]];
                    }
                }
            }
            catch(ArgumentOutOfRangeException)
            {
                throw new ArgumentOutOfRangeException("Sub matrix indices.");
            }
            return matrix;
        }

        /// <summary>
        /// Get a sub matrix
        /// </summary>
        /// <param name="i0">Start row index</param>
        /// <param name="i1">End row index</param>
        /// <param name="c">Array of column indices</param>
        /// <returns>_array(i0:i1, c(:))</returns>
        /// <exception>ArgumentOutOfRangeException : sub matrix indices.</exception>
        public Matrix GetMatrix(int i0, int i1, int[] c)
        {
            var m = i1 - i0 + 1;
            var n = c.Length;
            var matrix = new Matrix(m, n);
            var array = matrix.GetArray();
            try
            {
                for (int i = i0; i <= i1; ++i)
                {
                    for (int j = 0; j < n; ++j)
                    {
                        array[i-i0][j] = _array[i][c[j]];
                    }
                }
            }
            catch (ArgumentOutOfRangeException)
            {
                throw new ArgumentOutOfRangeException("Sub matrix indices.");
            }
            return matrix;
        }

        /// <summary>
        /// Get a sub matrix
        /// </summary>
        /// <param name="r">Array of row indices</param>
        /// <param name="j0">Start column index</param>
        /// <param name="j1">End column index</param>
        /// <returns>_array(r(:), j0:j1)</returns>
        /// <exception>ArgumentOutOfRangeException : sub matrix indices.</exception>
        public Matrix GetMatrix(int[] r, int j0, int j1)
        {
            var m = r.Length;
            var n = j1 - j0 + 1;
            var matrix = new Matrix(m, n);
            var array = matrix.GetArray();
            try
            {
                for (int i = 0; i < m; ++i)
                {
                    for (int j = j0; j <= j1; ++j)
                    {
                        array[i][j - j0] = _array[r[i]][j];
                    }
                }
            }
            catch (ArgumentOutOfRangeException)
            {
                throw new ArgumentOutOfRangeException("Sub matrix indices.");
            }
            return matrix;
        }

        /// <summary>
        /// Set a sub matrix
        /// </summary>
        /// <param name="i0">Start row index</param>
        /// <param name="i1">End row index</param>
        /// <param name="j0">Start column index</param>
        /// <param name="j1">End column index</param>
        /// <param name="matrix">_array(i0:i1, j0:j1)</param>
        /// <exception>ArgumentOutOfRangeException : sub matrix indices.</exception>
        public void SetMatrix(int i0, int i1, int j0, int j1, Matrix matrix)
        {
            try
            {
                for (int i = i0; i <= i1; ++i)
                {
                    for (int j = j0; j <= j1; ++j)
                    {
                        _array[i][j] = matrix.Get(i - i0, j - j0);
                    }
                }
            }
            catch(ArgumentOutOfRangeException)
            {
                throw new ArgumentOutOfRangeException("Sub matrix indices.");
            }
        }

        /// <summary>
        /// Set a sub matrix
        /// </summary>
        /// <param name="r">Array of row indices</param>
        /// <param name="c">Array of column indices</param>
        /// <param name="matrix">_array(r(:), c(:))</param>
        /// <exception>ArgumentOutOfRangeException : sub matrix indices.</exception>
        public void SetMatrix(int[] r, int[] c, Matrix matrix)
        {
            try
            {
                for (int i = 0, l1 = r.Length; i < l1; ++i)
                {
                    for (int j = 0, l2 = c.Length; j < l2; ++j)
                    {
                        _array[r[i]][c[j]] = matrix.Get(i, j);
                    }
                }
            }
            catch (ArgumentOutOfRangeException)
            {
                throw new ArgumentOutOfRangeException("Sub matrix indices.");
            }
        }

        /// <summary>
        /// Set a sub matrix
        /// </summary>
        /// <param name="r">Array of row indices</param>
        /// <param name="j0">Start column index</param>
        /// <param name="j1">End column index</param>
        /// <param name="matrix">_array(r(:), j0:j1)</param>
        /// <exception>ArgumentOutOfRangeException : sub matrix indices.</exception>
        public void SetMatrix(int[] r, int j0, int j1, Matrix matrix)
        {
            try
            {
                for (int i = 0, l1 = r.Length; i < l1; ++i)
                {
                    for (int j = j0; j <= j1; ++j)
                    {
                        _array[r[i]][j] = matrix.Get(i, j - j0);
                    }
                }
            }
            catch (ArgumentOutOfRangeException)
            {
                throw new ArgumentOutOfRangeException("Sub matrix indices.");
            }
        }

        /// <summary>
        /// Set a sub matrix
        /// </summary>
        /// <param name="i0">Start row index</param>
        /// <param name="i1">End row index</param>
        /// <param name="c">Array of column indices</param>
        /// <param name="matrix">_array(i0:i1, c(:))</param>
        /// <exception>ArgumentOutOfRangeException : sub matrix indices.</exception>
        public void SetMatrix(int i0, int i1, int[] c, Matrix matrix)
        {
            try
            {
                for (int i = i0; i <= i1; ++i)
                {
                    for (int j = 0, l2 = c.Length; j < l2; ++j)
                    {
                        _array[i][c[j]] = matrix.Get(i - i0, j);
                    }
                }
            }
            catch (ArgumentOutOfRangeException)
            {
                throw new ArgumentOutOfRangeException("Sub matrix indices.");
            }
        }

        /// <summary>
        /// Matrix transpose
        /// </summary>
        /// <returns>A'</returns>
        public Matrix Transpose()
        {
            var r = new Matrix(column, row);
            var array = r.GetArray();
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    array[j][i] = _array[i][j];
                }
            }
            return r;
        }

        /// <summary>
        /// One norm
        /// </summary>
        /// <returns>maximum column sum</returns>
        public double Norm1()
        {
            double f = 0;
            for (int j = 0; j < column; ++j)
            {
                double s = 0;
                for (int i = 0; i < row; ++i)
                {
                    s += Math.Abs(_array[i][j]);
                }

                f = Math.Max(f, s);
            }

            return f;
        }

        /// <summary>
        /// Two norm
        /// </summary>
        /// <returns>maximum singular value</returns>
        public double Norm2()
        {
            return (new SingularValueDecomposition(this).Norm2());
        }

        /// <summary>
        /// Infinity norm
        /// </summary>
        /// <returns>maximum row sum</returns>
        public double NormInf()
        {
            double f = 0;
            for (int i = 0; i < row; ++i)
            {
                double s = 0;
                for (int j = 0; j < column; ++j)
                {
                    s += Math.Abs(_array[i][j]);
                }

                f = Math.Max(f, s);
            }

            return f;
        }

        /// <summary>
        /// Frobenius norm
        /// </summary>
        /// <returns>sqrt of sum of squares of all elements</returns>
        public double NormF()
        {
            double f = 0;
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    f = Maths.Hypot(f, _array[i][j]);
                }
            }
            return f;
        }

        /// <summary>
        /// A = A + B
        /// </summary>
        /// <param name="b"></param>
        /// <returns>A += B</returns>
        public Matrix PlusEquals(Matrix b)
        {
            CheckMatrixDimensions(this, b);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    _array[i][j] += b._array[i][j];
                }
            }
            return this;
        }

        /// <summary>
        /// A = A - B
        /// </summary>
        /// <param name="b"></param>
        /// <returns>A -= B</returns>
        public Matrix MinusEquals(Matrix b)
        {
            CheckMatrixDimensions(this, b);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    _array[i][j] -= b._array[i][j];
                }
            }
            return this;
        }


        /// <summary>
        /// Element-by-element multiplication, C = A .* B
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>A.*B</returns>
        public Matrix ArrayTimes(Matrix b)
        {
            CheckMatrixDimensions(this, b);
            Matrix x = new Matrix(row, column);
            var array = x.GetArray();
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    array[i][j] = _array[i][j] * b._array[i][j];
                }
            }
            return x;
        }

        /// <summary>
        /// Element-by-element multiplication in place, A = A .* B
        /// </summary>
        /// <param name="b">another matrix</param>
        /// <returns>A.*B</returns>
        public Matrix ArrayTimesEquals(Matrix b)
        {
            CheckMatrixDimensions(this, b);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    _array[i][j] *= b._array[i][j];
                }
            }
            return this;
        }

        /// <summary>
        /// Element-by-element right division, C = A ./ B
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>A./B</returns>
        public Matrix ArrayRightDivide(Matrix b)
        {
            CheckMatrixDimensions(this, b);
            Matrix x = new Matrix(row, column);
            var array = x.GetArray();
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    array[i][j] = _array[i][j] / b._array[i][j];
                }
            }
            return x;
        }

        /// <summary>
        /// Element-by-element right division in place, A = A./B
        /// </summary>
        /// <param name="b">another matrix</param>
        /// <returns>A./B</returns>
        public Matrix ArrayRightDivideEquals(Matrix b)
        {
            CheckMatrixDimensions(this, b);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    _array[i][j] /= b._array[i][j];
                }
            }
            return this;
        }

        /// <summary>
        /// Element-by-element left division, C = B./A
        /// </summary>
        /// <param name="b"></param>
        /// <returns>B./A</returns>
        public Matrix ArrayLeftDivide(Matrix b)
        {
            return b.ArrayRightDivide(this);
        }

        /// <summary>
        /// Element-by-element left division in place, C = B./A
        /// </summary>
        /// <param name="b"></param>
        /// <returns>B./A</returns>
        public Matrix ArrayLeftDivideEquals(Matrix b)
        {
            CheckMatrixDimensions(this, b);
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    _array[i][j] = b._array[i][j] / _array[i][j];
                }
            }
            return this;
        }

        /// <summary>
        /// Multiply a matrix by a scalar in place, A = s*A
        /// </summary>
        /// <param name="s"></param>
        /// <returns>s*A</returns>
        public Matrix TimesEquals(double s)
        {
            for (int i = 0; i < row; ++i)
            {
                for (int j = 0; j < column; ++j)
                {
                    _array[i][j] *= s;
                }
            }

            return this;
        }

        /// <summary>
        /// Linear algebraic matrix multiplication, A * B
        /// </summary>
        /// <param name="b">another matrix</param>
        /// <returns>Matrix product, A * B</returns>
        /// <exception>ArgumentException : Matrix inner dimensions must agree.</exception>
        public Matrix Times(Matrix b)
        {
            if (b.row != column)
            {
                throw new ArgumentException("Matrix inner dimensions must agree.");
            }

            Matrix x = new Matrix(row, b.column);
            var array = x.GetArray();
            var bcolj = new double[column];
            for (int j = 0; j < b.column; ++j)
            {
                for (int k = 0; k < column; ++k)
                {
                    bcolj[k] = b._array[k][j];
                }
                for (int i = 0; i < row; ++i)
                {
                    var arowi = _array[i];
                    double sum = 0;
                    for (int k = 0; k < column; ++k)
                    {
                        sum += arowi[k] * bcolj[k];
                    }
                    array[i][j] = sum;
                }
            }

            return x;
        }

        /// <summary>
        /// LU Decomposition
        /// </summary>
        /// <return>LUDecomposition</return>
        public LUDecomposition LU()
        {
            return new LUDecomposition(this);
        }

        /// <summary>
        /// QR Decomposition
        /// </summary>
        /// <return>QRDecomposition</return>
        public QRDecomposition QR()
        {
            return new QRDecomposition(this);
        }

        /// <summary>
        /// Cholesky Decomposition
        /// </summary>
        /// <returns>CholeskyDecomposition</returns>
        public CholeskyDecomposition Cholesky()
        {
            return new CholeskyDecomposition(this);
        }

        /// <summary>
        /// Singular Value Decomposition
        /// </summary>
        /// <returns>SingularValueDecomposition</returns>
        public SingularValueDecomposition SVD()
        {
            return new SingularValueDecomposition(this);
        }

        /// <summary>
        /// Eigenvalue Decomposition
        /// </summary>
        /// <returns>EigenvalueDecomposition</returns>
        public EigenvalueDecomposition Eig()
        {
            return new EigenvalueDecomposition(this);
        }

        /// <summary>
        /// Solve A*X = B
        /// </summary>
        /// <param name="b">right hand side</param>
        /// <returns>solution if A is square, least squares solution otherwise</returns>
        public Matrix Solve(Matrix b)
        {
            return (row == column ? new LUDecomposition(this).Solve(b) :
                               new QRDecomposition(this).Solve(b));
        }

        /// <summary>
        /// Solve X*A = B, which is alse A'*X' = B'
        /// </summary>
        /// <param name="b">right hand side</param>
        /// <returns>solution if A is square, least squares solution otherwise.</returns>
        public Matrix SolveTranspose(Matrix b)
        {
            return Transpose().Solve(b.Transpose());
        }

        /// <summary>
        /// Matrix inverse or pseudo inverse
        /// </summary>
        /// <returns>inverse(A) if A is square, pseudoinverse otherwise.</returns>
        public Matrix Inverse()
        {
            return Solve(Identity(row, row));
        }

        /// <summary>
        /// Matrix determinant
        /// </summary>
        /// <returns>determinant</returns>
        public double Det()
        {
            return new LUDecomposition(this).Det();
        }

        /// <summary>
        /// Matrix rank
        /// </summary>
        /// <returns>effective numerical rank, obtained from SVD.</returns>
        public int Rank()
        {
            return new SingularValueDecomposition(this).Rank();
        }

        /// <summary>
        /// Matrix condition (2 norm)
        /// </summary>
        /// <returns>ratio of largest to smallest singular value.</returns>
        public double Cond()
        {
            return new SingularValueDecomposition(this).Cond();
        }

        /// <summary>
        /// Matrix trace
        /// </summary>
        /// <returns>sum of the diagonal elements.</returns>
        public double Trace()
        {
            double t = 0;
            var length = Math.Min(row, column);
            for (int i = 0; i < length; ++i)
            {
                t += _array[i][i];
            }

            return t;
        }

        #endregion

        #region operator
        /// <summary>
        /// Unary minus
        /// </summary>
        /// <param name="a"></param>
        /// <returns>-A</returns>
        public static Matrix operator -(Matrix a)
        {
            Matrix r = new Matrix(a.row, a.column);
            for (int i = 0; i < a.row; ++i)
            {
                for (int j = 0; j < a.column; ++j)
                {
                    r._array[i][j] = -a._array[i][j];
                }
            }

            return r;
        }

        /// <summary>
        /// C = A + B
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>A + B</returns>
        public static Matrix operator+(Matrix a, Matrix b)
        {
            CheckMatrixDimensions(a, b);
            Matrix r = new Matrix(a.row, a.column);
            for (int i = 0; i < a.row; ++i)
            {
                for (int j = 0; j < a.column; ++j)
                {
                    r._array[i][j] = a._array[i][j] + b._array[i][j];
                }
            }

            return r;
        }

        /// <summary>
        /// C = A - B
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <returns>A - B</returns>
        public static Matrix operator-(Matrix a, Matrix b)
        {
            CheckMatrixDimensions(a, b);
            Matrix r = new Matrix(a.row, a.column);
            for (int i = 0; i < a.row; ++i)
            {
                for (int j = 0; j < a.column; ++j)
                {
                    r._array[i][j] = a._array[i][j] - b._array[i][j];
                }
            }

            return r;
        }

        /// <summary>
        /// Multiply a matrix by a scalar, C = A*s
        /// </summary>
        /// <param name="a"></param>
        /// <param name="s"></param>
        /// <returns>A*s</returns>
        public static Matrix operator*(Matrix a, double s)
        {
            var x = new Matrix(a.row, a.column);
            var array = x.GetArray();
            for (int i = 0; i < a.row; ++i)
            {
                for (int j = 0; j < a.column; ++j)
                {
                    array[i][j] *= s;
                }
            }

            return x;
        }

        /// <summary>
        /// Multiply a matrix by a scalar, C = s*A
        /// </summary>
        /// <param name="a"></param>
        /// <param name="s"></param>
        /// <returns>s*A</returns>
        public static Matrix operator *(double s, Matrix a)
        {
            return a * s;
        }

        #endregion

        #region Private methods
        /// <summary>
        /// Check if size(A) == size(B)
        /// </summary>
        /// <param name="b"></param>
        private static void CheckMatrixDimensions(Matrix a, Matrix b)
        {
            if (b.row != a.row || a.column != b.column)
            {
                throw new ArgumentException("Matrix dimensions must agree.");
            }
        }
        #endregion

        #region Static methods
        /// <summary>
        /// Construct a matrix from a copy of a 2-D array.
        /// </summary>
        /// <param name="original">Two-dimensional array of T.</param>
        /// <exception>ArgumentException : All rows must have the same length.</exception>
        public static Matrix ConstructWithCopy(double[][] original)
        {
            int m = original.Length;
            int n = original[0].Length;
            var matrix = new Matrix(m, n);
            var array = matrix.GetArray();
            for (int i = 0; i < m; ++i)
            {
                if (original[i].Length != n)
                {
                    throw new ArgumentException("All rows must have the same length.");
                }

                for (int j = 0; j < n; ++j)
                {
                    array[i][j] = original[i][j];
                }
            }

            return matrix;
        }

        private static Random _random = new Random();
        /// <summary>
        /// Generate matrix with random elements
        /// </summary>
        /// <param name="m">Number of rows</param>
        /// <param name="n">Number of colums</param>
        /// <returns>An m-by-n matrix with uniformly distributed random elements.</returns>
        public static Matrix Random(int m, int n)
        {
            Matrix x = new Matrix(m, n);
            var array = x.GetArray();
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    array[i][j] = _random.NextDouble();
                }
            }
            return x;
        }

        /// <summary>
        /// Generate identity matrix
        /// </summary>
        /// <param name="m">Number of rows</param>
        /// <param name="n">Number of colums</param>
        /// <returns>An m-by-n matrix with ones on the diagonal and zeros elsewhere</returns>
        public static Matrix Identity(int m, int n)
        {
            var x = new Matrix(m, n);
            var array = x.GetArray();
            for (int i = 0; i < m; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    array[i][j] = i == j ? 1 : 0;
                }
            }

            return x;
        }
        #endregion
    }
}
