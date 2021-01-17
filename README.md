# Csma
It from JAMA(Java Matrix Package),https://math.nist.gov/javanumerics/jama/

Csma = CSharp Matrix class.
<P>
    The CSharp Matrix Class provides the fundamental operations of numerical
    linear algebra. Various constructors create Matrices from two dimensional
    arrays of double precision floating point numbers. Various "gets" and 
    "sets" provide access to sub matrices and matrix elements. Several methods
    implement basic matrix arithmetic, including matrix addition and
    multiplication, matrix norms, and element-by-element array operations.
    All the operations in this version of the Matrix Class involve real matrices.
    Complex matrices may be handled in a future version.
<P>
    Five fundamental matrix decompositions, which consist of pairs or triples
    of matrices, permutation vectors, and the like, produce results in five
    decomposition classes. These decompositions are accessed by the Matrix
    class to compute solutions of simultaneous linear equations, determinants,
    inverses and other matrix functions. The five decompositions are:
<P>
<UL>
   <LI>Cholesky Decomposition of symmetric, positive definite matrices.
   <LI>LU Decomposition of rectangular matrices.
   <LI>QR Decomposition of rectangular matrices.
   <LI>Singular Value Decomposition of retangular matrices.
   <LI>Eigenvalue Decomposition of both symmetric and nonsymmetric square matrices.
</UL>
<DL>
<DT><B>Example of use:</B></DT>
<P>
<DD>Solve a linear system A x = b and compute the residual norm,||b - A x||.
<P><PRE>
      double[][] vals = {{1,2,3},{4,5,6},{7,8,9}};
      Matrix A = new Matrix(vals);
      Matrix b = Matrix.Random(3, 1);
      Matrix x = A.Solve(b);
      Matrix r = A.Times(x) - b;
      double rnorm = r.NormInf();
</PRE></DD>
</DL>
