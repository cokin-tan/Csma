using System;

namespace Csma
{
    public static class Test
    {
        public static void StartTest()
        {
            Matrix A, B, C, Z, O, I, R, S, X, SUB, M, T, SQ, DEF, SOL;
            // Uncomment this to test IO in a different locale.
            // Locale.setDefault(Locale.GERMAN);
            int errorCount = 0;
            int warningCount = 0;
            double tmp, s;
            double[] columnwise = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12 };
            double[] rowwise = { 1, 4, 7, 10, 2, 5, 8, 11, 3, 6, 9, 12 };
            double[][] avals = { new double[] { 1, 4, 7, 10 }, new double[] { 2, 5, 8, 11 }, new double[] { 3, 6, 9, 12 } };
            double[][] rankdef = avals;
            double[][] tvals = { new double[] { 1, 2, 3 }, new double[] { 4, 5, 6 }, new double[] { 7, 8, 9 }, new double[] { 10, 11, 12 } };
            double[][] subavals = { new double[] { 5, 8, 11 }, new double[] { 6, 9, 12 } };
            double[][] rvals = { new double[] { 1, 4, 7 }, new double[] { 2, 5, 8, 11 }, new double[] { 3, 6, 9, 12 } };
            double[][] pvals = { new double[] { 4, 1, 1 }, new double[] { 1, 2, 3 }, new double[] { 1, 3, 6 } };
            double[][] ivals = { new double[] { 1, 0, 0, 0 }, new double[] { 0, 1, 0, 0 }, new double[] { 0, 0, 1, 0 } };
            double[][] evals = { new double[] { 0, 1, 0, 0 }, new double[] { 1, 0, 2e-7, 0 }, new double[] { 0, -2e-7, 0, 1 }, new double[] { 0, 0, 1, 0 } };
            double[][] square = { new double[] { 166, 188, 210 }, new double[] { 188, 214, 240 }, new double[] { 210, 240, 270 } };
            double[][] sqSolution = { new double[] { 13 }, new double[] { 15 } };
            double[][] condmat = { new double[] { 1, 3 }, new double[] { 7, 9 } };
            int rows = 3, cols = 4;
            int invalidld = 5;/* should trigger bad shape for construction with val */
            int raggedr = 0; /* (raggedr,raggedc) should be out of bounds in ragged array */
            int raggedc = 4;
            int validld = 3; /* leading dimension of intended test Matrices */
            int nonconformld = 4; /* leading dimension which is valid, but nonconforming */
            int ib = 1, ie = 2, jb = 1, je = 3; /* index ranges for sub Matrix */
            int[] rowindexset = { 1, 2 };
            int[] badrowindexset = { 1, 3 };
            int[] columnindexset = { 1, 2, 3 };
            int[] badcolumnindexset = { 1, 2, 4 };
            double columnsummax = 33;
            double rowsummax = 30;
            double sumofdiagonals = 15;
            double sumofsquares = 650;

            /** 
             * Constructors and constructor-like methods:
             * double[], int
             * double[][]
             * int, int
             * int, int, double
             * int, int, double[][]
             * constructWithCopy(double[][])
             * random(int,int)
             * identity(int)
             * **/
            Console.WriteLine("\nTesting constructors and constructor-like methods...\n");
            try
            {
                A = new Matrix(columnwise, invalidld);
                errorCount = try_failure(errorCount, "Catch invalid length in packed constructor... ",
                     "exception not thrown for invalid input");
            }
            catch (Exception e)
            {
                try_success("Catch invalid length in packed constructor... ",
                     e.Message);
            }

            try
            {
                /** check that exception is thrown in default constructor
                    if input array is 'ragged' **/
                A = new Matrix(rvals);
                tmp = A.Get(raggedr, raggedc);
            }
            catch (ArgumentException e)
            {
                try_success("Catch ragged input to default constructor... ",
                             e.Message);
            }
            catch (IndexOutOfRangeException e)
            {
                errorCount = try_failure(errorCount, "Catch ragged input to constructor... ",
                            "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
            }

            try
            {
                /** check that exception is thrown in constructWithCopy
                    if input array is 'ragged' **/
                A = Matrix.ConstructWithCopy(rvals);
                tmp = A.Get(raggedr, raggedc);
            }
            catch (ArgumentException e)
            {
                try_success("Catch ragged input to constructWithCopy... ", e.Message);
            }
            catch (IndexOutOfRangeException e)
            {
                errorCount = try_failure(errorCount, "Catch ragged input to constructWithCopy... ", "exception not thrown in construction...ArrayIndexOutOfBoundsException thrown later");
            }

            A = new Matrix(columnwise, validld);
            B = new Matrix(avals);
            tmp = B.Get(0, 0);
            avals[0][0] = 0.0;
            C = B - A;
            avals[0][0] = tmp;
            B = Matrix.ConstructWithCopy(avals);
            tmp = B.Get(0, 0);
            avals[0][0] = 0.0;
            if ((tmp - B.Get(0, 0)) != 0.0)
            {
                /** check that constructWithCopy behaves properly **/
                errorCount = try_failure(errorCount, "constructWithCopy... ", "copy not effected... data visible outside");
            }
            else
            {
                try_success("constructWithCopy... ", "");
            }
            avals[0][0] = columnwise[0];
            I = new Matrix(ivals);
            try
            {
                check(I, Matrix.Identity(3, 4));
                try_success("identity... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "identity... ", "identity Matrix not successfully created");
            }
            /**   
             * Access Methods:
             * getColumnDimension()
             * getRowDimension()
             * getArray()
             * getArrayCopy()
             * getColumnPackedCopy()
             * getRowPackedCopy()
             * get(int,int)
             * getMatrix(int,int,int,int)
             * getMatrix(int,int,int[])
             * getMatrix(int[],int,int)
             * getMatrix(int[],int[])
             * set(int,int,double)
             * setMatrix(int,int,int,int,Matrix)
             * setMatrix(int,int,int[],Matrix)
             * setMatrix(int[],int,int,Matrix)
             * setMatrix(int[],int[],Matrix)
             * **/
            Console.WriteLine("\nTesting access methods...\n");
            /**
             * Various get methods:
             * **/
            B = new Matrix(avals);
            if (B.row != rows)
            {
                errorCount = try_failure(errorCount, "getRowDimension... ", "");
            }
            else
            {
                try_success("getRowDimension... ", "");
            }
            if (B.column != cols)
            {
                errorCount = try_failure(errorCount, "getColumnDimension... ", "");
            }
            else
            {
                try_success("getColumnDimension... ", "");
            }
            B = new Matrix(avals);
            double[][] barray = B.GetArray();
            if (barray != avals)
            {
                errorCount = try_failure(errorCount, "getArray... ", "");
            }
            else
            {
                try_success("getArray... ", "");
            }
            barray = B.GetArrayCopy();
            if (barray == avals)
            {
                errorCount = try_failure(errorCount, "getArrayCopy... ", "data not (deep) copied");
            }
            try
            {
                check(barray, avals);
                try_success("getArrayCopy... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "getArrayCopy... ", "data not successfully (deep) copied");
            }

            double[] bpacked = B.GetColumnPackedCopy();
            try
            {
                check(bpacked, columnwise);
                try_success("getColumnPackedCopy... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "getColumnPackedCopy... ", "data not successfully (deep) copied by columns");
            }
            bpacked = B.GetRowPackedCopy();
            try
            {
                check(bpacked, rowwise);
                try_success("getRowPackedCopy... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "getRowPackedCopy... ", "data not successfully (deep) copied by rows");
            }

            try
            {
                tmp = B.Get(B.row, B.column - 1);
                errorCount = try_failure(errorCount, "get(int,int)... ", "OutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    tmp = B.Get(B.row - 1, B.column);
                    errorCount = try_failure(errorCount, "get(int,int)... ", "OutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("get(int,int)... OutofBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "get(int,int)... ", "OutOfBoundsException expected but not thrown");
            }

            try
            {
                if (B.Get(B.row - 1, B.column - 1) !=
                    avals[B.row - 1][B.column - 1])
                {
                    errorCount = try_failure(errorCount, "get(int,int)... ", "Matrix entry (i,j) not successfully retreived");
                }
                else
                {
                    try_success("get(int,int)... ", "");
                }
            }
            catch (IndexOutOfRangeException e)
            {
                errorCount = try_failure(errorCount, "get(int,int)... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            SUB = new Matrix(subavals);
            try
            {
                M = B.GetMatrix(ib, ie + B.row + 1, jb, je);
                errorCount = try_failure(errorCount, "getMatrix(int,int,int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    M = B.GetMatrix(ib, ie, jb, je + B.column + 1);
                    errorCount = try_failure(errorCount, "getMatrix(int,int,int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("getMatrix(int,int,int,int)... ArrayIndexOutOfBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "getMatrix(int,int,int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }

            try
            {
                M = B.GetMatrix(ib, ie, jb, je);
                try
                {
                    check(SUB, M);
                    try_success("getMatrix(int,int,int,int)... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "getMatrix(int,int,int,int)... ", "submatrix not successfully retreived");
                }
            }
            catch (IndexOutOfRangeException e)
            {
                errorCount = try_failure(errorCount, "getMatrix(int,int,int,int)... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            try
            {
                M = B.GetMatrix(ib, ie, badcolumnindexset);
                errorCount = try_failure(errorCount, "getMatrix(int,int,int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    M = B.GetMatrix(ib, ie + B.row + 1, columnindexset);
                    errorCount = try_failure(errorCount, "getMatrix(int,int,int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("getMatrix(int,int,int[])... ArrayIndexOutOfBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "getMatrix(int,int,int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }

            try
            {
                M = B.GetMatrix(ib, ie, columnindexset);
                try
                {
                    check(SUB, M);
                    try_success("getMatrix(int,int,int[])... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "getMatrix(int,int,int[])... ", "submatrix not successfully retreived");
                }
            }
            catch (IndexOutOfRangeException e)
            {
                errorCount = try_failure(errorCount, "getMatrix(int,int,int[])... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            try
            {
                M = B.GetMatrix(badrowindexset, jb, je);
                errorCount = try_failure(errorCount, "getMatrix(int[],int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    M = B.GetMatrix(rowindexset, jb, je + B.column + 1);
                    errorCount = try_failure(errorCount, "getMatrix(int[],int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("getMatrix(int[],int,int)... ArrayIndexOutOfBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "getMatrix(int[],int,int)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }

            try
            {
                M = B.GetMatrix(rowindexset, jb, je);
                try
                {
                    check(SUB, M);
                    try_success("getMatrix(int[],int,int)... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "getMatrix(int[],int,int)... ", "submatrix not successfully retreived");
                }
            }
            catch (IndexOutOfRangeException e)
            {
                errorCount = try_failure(errorCount, "getMatrix(int[],int,int)... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            try
            {
                M = B.GetMatrix(badrowindexset, columnindexset);
                errorCount = try_failure(errorCount, "getMatrix(int[],int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    M = B.GetMatrix(rowindexset, badcolumnindexset);
                    errorCount = try_failure(errorCount, "getMatrix(int[],int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("getMatrix(int[],int[])... ArrayIndexOutOfBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "getMatrix(int[],int[])... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }

            try
            {
                M = B.GetMatrix(rowindexset, columnindexset);
                try
                {
                    check(SUB, M);
                    try_success("getMatrix(int[],int[])... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "getMatrix(int[],int[])... ", "submatrix not successfully retreived");
                }
            }
            catch (IndexOutOfRangeException e)
            {
                errorCount = try_failure(errorCount, "getMatrix(int[],int[])... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            /**
             * Various set methods:
             * **/
            try
            {
                B.Set(B.row, B.column - 1, 0);
                errorCount = try_failure(errorCount, "set(int,int,double)... ", "OutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    B.Set(B.row - 1, B.column, 0);
                    errorCount = try_failure(errorCount, "set(int,int,double)... ", "OutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("set(int,int,double)... OutofBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "set(int,int,double)... ", "OutOfBoundsException expected but not thrown");
            }

            try
            {
                B.Set(ib, jb, 0);
                tmp = B.Get(ib, jb);
                try
                {
                    check(tmp, 0);
                    try_success("set(int,int,double)... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "set(int,int,double)... ", "Matrix element not successfully set");
                }
            }
            catch (IndexOutOfRangeException e1)
            {
                errorCount = try_failure(errorCount, "set(int,int,double)... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            M = new Matrix(2, 3, 0);
            try
            {
                B.SetMatrix(ib, ie + B.row + 1, jb, je, M);
                errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    B.SetMatrix(ib, ie, jb, je + B.column + 1, M);
                    errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("setMatrix(int,int,int,int,Matrix)... ArrayIndexOutOfBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }

            try
            {
                B.SetMatrix(ib, ie, jb, je, M);
                try
                {
                    check(M - (B.GetMatrix(ib, ie, jb, je)), M);
                    try_success("setMatrix(int,int,int,int,Matrix)... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ", "submatrix not successfully set");
                }
                B.SetMatrix(ib, ie, jb, je, SUB);
            }
            catch (IndexOutOfRangeException e1)
            {
                errorCount = try_failure(errorCount, "setMatrix(int,int,int,int,Matrix)... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            try
            {
                B.SetMatrix(ib, ie + B.row + 1, columnindexset, M);
                errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    B.SetMatrix(ib, ie, badcolumnindexset, M);
                    errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("setMatrix(int,int,int[],Matrix)... ArrayIndexOutOfBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }

            try
            {
                B.SetMatrix(ib, ie, columnindexset, M);
                try
                {
                    check(M - (B.GetMatrix(ib, ie, columnindexset)), M);
                    try_success("setMatrix(int,int,int[],Matrix)... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ", "submatrix not successfully set");
                }
                B.SetMatrix(ib, ie, jb, je, SUB);
            }
            catch (IndexOutOfRangeException e1)
            {
                errorCount = try_failure(errorCount, "setMatrix(int,int,int[],Matrix)... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            try
            {
                B.SetMatrix(rowindexset, jb, je + B.column + 1, M);
                errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    B.SetMatrix(badrowindexset, jb, je, M);
                    errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("setMatrix(int[],int,int,Matrix)... ArrayIndexOutOfBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }

            try
            {
                B.SetMatrix(rowindexset, jb, je, M);
                try
                {
                    check(M - (B.GetMatrix(rowindexset, jb, je)), M);
                    try_success("setMatrix(int[],int,int,Matrix)... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ", "submatrix not successfully set");
                }
                B.SetMatrix(ib, ie, jb, je, SUB);
            }
            catch (IndexOutOfRangeException e1)
            {
                errorCount = try_failure(errorCount, "setMatrix(int[],int,int,Matrix)... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            try
            {
                B.SetMatrix(rowindexset, badcolumnindexset, M);
                errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }
            catch (IndexOutOfRangeException e)
            {
                try
                {
                    B.SetMatrix(badrowindexset, columnindexset, M);
                    errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
                }
                catch (IndexOutOfRangeException e1)
                {
                    try_success("setMatrix(int[],int[],Matrix)... ArrayIndexOutOfBoundsException... ", "");
                }
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ", "ArrayIndexOutOfBoundsException expected but not thrown");
            }

            try
            {
                B.SetMatrix(rowindexset, columnindexset, M);
                try
                {
                    check(M - (B.GetMatrix(rowindexset, columnindexset)), M);
                    try_success("setMatrix(int[],int[],Matrix)... ", "");
                }
                catch (Exception e)
                {
                    errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ", "submatrix not successfully set");
                }
            }
            catch (IndexOutOfRangeException e1)
            {
                errorCount = try_failure(errorCount, "setMatrix(int[],int[],Matrix)... ", "Unexpected ArrayIndexOutOfBoundsException");
            }

            /** 
             * Array-like methods:
             * minus
             * minusEquals
             * plus
             * plusEquals
             * arrayLeftDivide
             * arrayLeftDivideEquals
             * arrayRightDivide
             * arrayRightDivideEquals
             * arrayTimes
             * arrayTimesEquals
             * uminus
             * **/
            Console.WriteLine("\nTesting array-like methods...\n");

            S = new Matrix(columnwise, nonconformld);
            R = Matrix.Random(A.row, A.column);
            A = R;
            try
            {
                S = A - (S);
                errorCount = try_failure(errorCount, "minus conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("minus conformance check... ", "");
            }

            if ((A - (R)).Norm1() != 0)
            {
                errorCount = try_failure(errorCount, "minus... ", "(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)");
            }
            else
            {
                try_success("minus... ", "");
            }

            A = R.Copy();
            A.MinusEquals(R);
            Z = new Matrix(A.row, A.column);
            try
            {
                A.MinusEquals(S);
                errorCount = try_failure(errorCount, "minusEquals conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("minusEquals conformance check... ", "");
            }

            if ((A - (Z)).Norm1() != 0)
            {
                errorCount = try_failure(errorCount, "minusEquals... ", "(difference of identical Matrices is nonzero,\nSubsequent use of minus should be suspect)");
            }
            else
            {
                try_success("minusEquals... ", "");
            }

            A = R.Copy();
            B = Matrix.Random(A.row, A.column);
            C = A - (B);
            try
            {
                S = A + (S);
                errorCount = try_failure(errorCount, "plus conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("plus conformance check... ", "");
            }

            try
            {
                check(C + (B), A);
                try_success("plus... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "plus... ", "(C = A - B, but C + B != A)");
            }

            C = A - (B);
            C.PlusEquals(B);
            try
            {
                A.PlusEquals(S);
                errorCount = try_failure(errorCount, "plusEquals conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("plusEquals conformance check... ", "");
            }

            try
            {
                check(C, A);
                try_success("plusEquals... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "plusEquals... ", "(C = A - B, but C = C + B != A)");
            }

            A = -R;
            try
            {
                check(A + (R), Z);
                try_success("uminus... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "uminus... ", "(-A + A != zeros)");
            }

            A = R.Copy();
            O = new Matrix(A.row, A.column, 1.0);
            C = A.ArrayLeftDivide(R);
            try
            {
                S = A.ArrayLeftDivide(S);
                errorCount = try_failure(errorCount, "arrayLeftDivide conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("arrayLeftDivide conformance check... ", "");
            }

            try
            {
                check(C, O);
                try_success("arrayLeftDivide... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "arrayLeftDivide... ", "(M.\\M != ones)");
            }

            try
            {
                A.ArrayLeftDivideEquals(S);
                errorCount = try_failure(errorCount, "arrayLeftDivideEquals conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("arrayLeftDivideEquals conformance check... ", "");
            }

            A.ArrayLeftDivideEquals(R);
            try
            {
                check(A, O);
                try_success("arrayLeftDivideEquals... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "arrayLeftDivideEquals... ", "(M.\\M != ones)");
            }

            A = R.Copy();
            try
            {
                A.ArrayRightDivide(S);
                errorCount = try_failure(errorCount, "arrayRightDivide conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("arrayRightDivide conformance check... ", "");
            }

            C = A.ArrayRightDivide(R);
            try
            {
                check(C, O);
                try_success("arrayRightDivide... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "arrayRightDivide... ", "(M./M != ones)");
            }

            try
            {
                A.ArrayRightDivideEquals(S);
                errorCount = try_failure(errorCount, "arrayRightDivideEquals conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("arrayRightDivideEquals conformance check... ", "");
            }

            A.ArrayRightDivideEquals(R);
            try
            {
                check(A, O);
                try_success("arrayRightDivideEquals... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "arrayRightDivideEquals... ", "(M./M != ones)");
            }

            A = R.Copy();
            B = Matrix.Random(A.row, A.column);
            try
            {
                S = A.ArrayTimes(S);
                errorCount = try_failure(errorCount, "arrayTimes conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("arrayTimes conformance check... ", "");
            }

            C = A.ArrayTimes(B);
            try
            {
                check(C.ArrayRightDivideEquals(B), A);
                try_success("arrayTimes... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "arrayTimes... ", "(A = R, C = A.*B, but C./B != A)");
            }

            try
            {
                A.ArrayTimesEquals(S);
                errorCount = try_failure(errorCount, "arrayTimesEquals conformance check... ", "nonconformance not raised");
            }
            catch (ArgumentException e)
            {
                try_success("arrayTimesEquals conformance check... ", "");
            }
            A.ArrayTimesEquals(B);
            try
            {
                check(A.ArrayRightDivideEquals(B), R);
                try_success("arrayTimesEquals... ", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "arrayTimesEquals... ", "(A = R, A = A.*B, but A./B != R)");
            }

            /** 
             * LA methods:
             * transpose
             * times
             * cond
             * rank
             * det
             * trace
             * norm1
             * norm2
             * normF
             * normInf
             * solve
             * solveTranspose
             * inverse
             * chol
             * eig
             * lu
             * qr
             * svd 
             * **/
            Console.WriteLine("\nTesting linear algebra methods...\n");
            A = new Matrix(columnwise, 3);
            T = new Matrix(tvals);
            T = A.Transpose();
            try
            {
                check(A.Transpose(), T);
                try_success("transpose...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "transpose()...", "transpose unsuccessful");
            }

            A.Transpose();
            try
            {
                check(A.Norm1(), columnsummax);
                try_success("norm1...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "norm1()...", "incorrect norm calculation");
            }

            try
            {
                check(A.NormInf(), rowsummax);
                try_success("normInf()...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "normInf()...", "incorrect norm calculation");
            }

            try
            {
                check(A.NormF(), Math.Sqrt(sumofsquares));
                try_success("normF...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "normF()...", "incorrect norm calculation");
            }

            try
            {
                check(A.Trace(), sumofdiagonals);
                try_success("trace()...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "trace()...", "incorrect trace calculation");
            }

            try
            {
                check(A.GetMatrix(0, A.row - 1, 0, A.row - 1).Det(), 0);
                try_success("det()...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "det()...", "incorrect determinant calculation");
            }
            SQ = new Matrix(square);
            try
            {
                check(A.Times(A.Transpose()), SQ);
                try_success("times(Matrix)...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "times(Matrix)...", "incorrect Matrix-Matrix product calculation");
            }

            try
            {
                check(A * (0), Z);
                try_success("times(double)...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "times(double)...", "incorrect Matrix-scalar product calculation");
            }

            A = new Matrix(columnwise, 4);
            QRDecomposition QR = A.QR();
            R = QR.GetR();
            try
            {
                check(A, QR.GetQ().Times(R));
                try_success("QRDecomposition...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "QRDecomposition...", "incorrect QR decomposition calculation");
            }

            SingularValueDecomposition SVD = A.SVD();
            try
            {
                check(A, SVD.GetU().Times(SVD.GetS().Times(SVD.GetV().Transpose())));
                try_success("SingularValueDecomposition...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "SingularValueDecomposition...", "incorrect singular value decomposition calculation");
            }

            DEF = new Matrix(rankdef);
            try
            {
                check(DEF.Rank(), Math.Min(DEF.row, DEF.column) - 1);
                try_success("rank()...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "rank()...", "incorrect rank calculation");
            }

            B = new Matrix(condmat);
            SVD = B.SVD();
            double[] singularvalues = SVD.GetSingularValues();
            try
            {
                check(B.Cond(), singularvalues[0] / singularvalues[Math.Min(B.row, B.column) - 1]);
                try_success("cond()...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "cond()...", "incorrect condition number calculation");
            }

            int n = A.column;
            A = A.GetMatrix(0, n - 1, 0, n - 1);
            A.Set(0, 0, 0);
            LUDecomposition LU = A.LU();
            try
            {
                check(A.GetMatrix(LU.GetPivot(), 0, n - 1), LU.GetL().Times(LU.GetU()));
                try_success("LUDecomposition...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "LUDecomposition...", "incorrect LU decomposition calculation");
            }

            X = A.Inverse();
            try
            {
                check(A.Times(X), Matrix.Identity(3, 3));
                try_success("inverse()...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "inverse()...", "incorrect inverse calculation");
            }

            O = new Matrix(SUB.row, 1, 1.0);
            SOL = new Matrix(sqSolution);
            SQ = SUB.GetMatrix(0, SUB.row - 1, 0, SUB.row - 1);
            try
            {
                check(SQ.Solve(SOL), O);
                try_success("solve()...", "");
            }
            catch (ArgumentException e1)
            {
                errorCount = try_failure(errorCount, "solve()...", e1.Message);
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "solve()...", e.Message);
            }

            A = new Matrix(pvals);
            CholeskyDecomposition Chol = A.Cholesky();
            Matrix L = Chol.GetL();
            try
            {
                check(A, L.Times(L.Transpose()));
                try_success("CholeskyDecomposition...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "CholeskyDecomposition...", "incorrect Cholesky decomposition calculation");
            }

            X = Chol.Solve(Matrix.Identity(3, 3));
            try
            {
                check(A.Times(X), Matrix.Identity(3, 3));
                try_success("CholeskyDecomposition solve()...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "CholeskyDecomposition solve()...", "incorrect Choleskydecomposition solve calculation");
            }

            EigenvalueDecomposition Eig = A.Eig();
            Matrix D = Eig.GetD();
            Matrix V = Eig.GetV();
            try
            {
                check(A.Times(V), V.Times(D));
                try_success("EigenvalueDecomposition (symmetric)...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "EigenvalueDecomposition (symmetric)...", "incorrect symmetric Eigenvalue decomposition calculation");
            }

            A = new Matrix(evals);
            Eig = A.Eig();
            D = Eig.GetD();
            V = Eig.GetV();
            try
            {
                check(A.Times(V), V.Times(D));
                try_success("EigenvalueDecomposition (nonsymmetric)...", "");
            }
            catch (Exception e)
            {
                errorCount = try_failure(errorCount, "EigenvalueDecomposition (nonsymmetric)...", "incorrect nonsymmetric Eigenvalue decomposition calculation");
            }
        }

        /** Check magnitude of difference of scalars. **/

        private static void check(double x, double y)
        {
            double eps = Math.Pow(2.0, -52.0);
            if (x == 0 & Math.Abs(y) < 10 * eps) return;
            if (y == 0 & Math.Abs(x) < 10 * eps) return;
            if (Math.Abs(x - y) > 10 * eps * Math.Max(Math.Abs(x), Math.Abs(y)))
            {
                throw new Exception("The difference x-y is too large: x = " + x.ToString("f6") + "  y = " + y.ToString("f6"));
            }
        }

        /** Check norm of difference of "vectors". **/

        private static void check(double[] x, double[] y)
        {
            if (x.Length == y.Length)
            {
                for (int i = 0; i < x.Length; i++)
                {
                    check(x[i], y[i]);
                }
            }
            else
            {
                throw new Exception("Attempt to compare vectors of different lengths");
            }
        }

        /** Check norm of difference of arrays. **/

        private static void check(double[][] x, double[][] y)
        {
            Matrix A = new Matrix(x);
            Matrix B = new Matrix(y);
            check(A, B);
        }

        /** Check norm of difference of Matrices. **/

        private static void check(Matrix X, Matrix Y)
        {
            double eps = Math.Pow(2.0, -52.0);
            if ((X.Norm1() == 0) & (Y.Norm1() < 10 * eps)) return;
            if ((Y.Norm1() == 0) & (X.Norm1() < 10 * eps)) return;
            if ((X - Y).Norm1() > 1000 * eps * Math.Max(X.Norm1(), Y.Norm1()))
            {
                throw new Exception("The norm of (X-Y) is too large: " + (X - Y).Norm1().ToString("f6"));
            }
        }

        private static void try_success(string s, string e)
        {
            Console.Write(">    " + s + "success\n");
            if ("" != e)
            {
                Console.Write(">      Message: " + e + "\n");
            }
        }

        private static int try_failure(int count, String s, String e)
        {
            Console.Write(">    " + s + "*** failure ***\n>      Message: " + e + "\n");
            return ++count;
        }

        private static int try_warning(int count, String s, String e)
        {
            Console.Write(">    " + s + "*** warning ***\n>      Message: " + e + "\n");
            return ++count;
        }
    }
}
