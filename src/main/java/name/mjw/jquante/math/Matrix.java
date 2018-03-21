package name.mjw.jquante.math;

import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * A general NxM real matrix.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class Matrix implements Cloneable {

	/**
	 * Holds value of property matrix.
	 */
	protected double[][] data;

	/**
	 * Holds value of property rowCount.
	 */
	protected int rowCount;

	/**
	 * Holds value of property columnCount.
	 */
	protected int columnCount;

	/**
	 * Creates a new instance of NxM Matrix
	 * 
	 * @param n
	 *            the first dimension
	 * @param m
	 *            the second dimension
	 */
	public Matrix(int n, int m) {
		data = new double[n][m];

		rowCount = n;
		columnCount = m;
	}

	private Matrix(int n, int m, boolean allocateStorage) {
		if (allocateStorage)
			data = new double[n][m];

		rowCount = n;
		columnCount = m;
	}

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public Matrix(int n) {
		this(n, n);
	}

	/**
	 * Creates a new instance of Matrix, based on already allocated 2D array
	 * 
	 * @param a
	 *            the 2D array
	 */
	public Matrix(double[][] a) {
		data = a;

		rowCount = a.length;
		columnCount = a[0].length;
	}

	/**
	 * Creates a diagonal matrix from a Vector
	 * 
	 * @param vec
	 *            the Vector
	 */
	public Matrix(Vector vec) {
		this(vec.getSize());

		for (int i = 0; i < data.length; i++)
			data[i][i] = vec.getVector(i);
	}

	/**
	 * Simple matrix addition of two square matrices: this + b
	 * 
	 * @param b
	 *            the matrix to which to add
	 * @return the result Cij = (Aij + Bij)
	 */
	public Matrix add(Matrix b) {
		Matrix c = new Matrix(rowCount);

		for (int i = 0; i < rowCount; i++) {
			for (int j = 0; j < columnCount; j++) {
				c.data[i][j] = data[i][j] + b.data[i][j];
			}
		}

		return c;
	}

	/**
	 * Simple matrix subtraction of two square matrices: this - b
	 * 
	 * @param b
	 *            the matrix to which to subtract
	 * @return the result Cij = (Aij - Bij)
	 */
	public Matrix sub(Matrix b) {
		Matrix c = new Matrix(rowCount);

		for (int i = 0; i < rowCount; i++) {
			for (int j = 0; j < columnCount; j++) {
				c.data[i][j] = data[i][j] - b.data[i][j];
			}
		}

		return c;
	}

	/**
	 * Simple matrix multiplication of two matrices: this * b
	 * 
	 * @param b
	 *            the matrix to which to multiply
	 * @return the result Cij = Sum(Aik*Bkj)
	 */
	public Matrix multiply(Matrix b) {
		Matrix c = new Matrix(rowCount, b.columnCount, false);

		c.data = Arrays.stream(data)
				.map(r -> IntStream.range(0, b.data[0].length).parallel()
						.mapToDouble(
								i -> IntStream.range(0, b.data.length).mapToDouble(j -> r[j] * b.data[j][i]).sum())
						.toArray())
				.toArray(double[][]::new);

		return c;
	}

	/**
	 * Defines a vector '.' operation on two matrices, 'as if' flattened out as NxM
	 * size vector. Both the matrices need to of same size to give correct result.
	 * No checks are made to make sure this is the case.
	 * 
	 * @param b
	 *            the matrix to do a dot product
	 * @return the value of dot product
	 */
	public double dotProduct(Matrix b) {
		double res = 0.0;

		for (int i = 0; i < rowCount; i++) {
			for (int j = 0; j < columnCount; j++) {
				res += data[i][j] * b.data[i][j];
			}
		}

		return res;
	}

	/**
	 * Transpose a real matrix A (A')
	 * 
	 * @return the transposed instance
	 */
	public Matrix transpose() {
		Matrix matrixT = new Matrix(columnCount, rowCount);

		for (int i = 0; i < rowCount; i++) {
			for (int j = 0; j < columnCount; j++) {
				matrixT.data[j][i] = data[i][j];
			}
		}

		return matrixT;
	}

	/**
	 * Tr(this) = sum(diagonal(this)); Only valid if this a square matrix.
	 * 
	 * @return the trace of this matrix
	 */
	public double getTrace() {
		double tr = 0.0;

		for (int i = 0; i < rowCount; i++)
			tr += data[i][i];

		return tr;
	}

	/**
	 * make the current matrix an identity matrix, all diagonals as 1.0 and
	 * non-diagonals zero<br>
	 * only sensible if a square matrix
	 */
	public void makeIdentity() {
		for (int i = 0; i < rowCount; i++) {
			for (int j = 0; j < i; j++) {
				data[i][j] = data[j][i] = 0.0;
			}
			data[i][i] = 1.0;
		}
	}

	/**
	 * rootMeanSquare() - method to compute the root mean square of elements of a
	 * diagonal. This is one type of matrix norm.
	 * <code>sqrt(sum(A<sub>ii</sub>))</code>.
	 * 
	 * @return double - the root mean square of the elements of the diagonal.
	 */
	public double rootMeanSquare() {
		double sum = 0.0;

		for (int i = 0; i < rowCount; i++) {
			sum += data[i][i] * data[i][i]; // sum of squares
		}

		return Math.sqrt(sum / rowCount);
	}

	/**
	 * Absolute sum of off-diagonal elements.
	 * 
	 * @return Sum(a[i][j], i!=j)
	 */
	public double sumOffDiagonal() {
		double sum = 0.0;
		int i;
		int j;

		for (i = 0; i < rowCount - 1; i++) {
			for (j = i + 1; j < rowCount; j++) {
				sum += Math.abs(data[i][j]);
			}
		}

		return sum;
	}

	/**
	 * Getter for property matrix.
	 * 
	 * @return Value of property matrix.
	 */
	public double[][] getData() {
		return this.data;
	}

	/**
	 * Getter for property matrix.
	 * 
	 * @param i
	 *            index i
	 * @param j
	 *            index j
	 * 
	 * @return Value of property matrix.
	 */
	public double getMatrixAt(int i, int j) {
		return this.data[i][j];
	}

	/**
	 * Setter for property matrix.
	 * 
	 * @param matrix
	 *            New value of property matrix.
	 */
	public void setMatrix(double[][] matrix) {
		this.data = matrix;

		rowCount = matrix.length;
		columnCount = matrix[0].length;
	}

	/**
	 * Getter for property matrix.
	 * 
	 * @param i
	 *            index i
	 * @param j
	 *            index j
	 * @param value
	 *            the new value
	 */
	public void setMatrixAt(int i, int j, double value) {
		this.data[i][j] = value;
	}

	/**
	 * Convert this Matrix object to a Vector object
	 * 
	 * @return the "flattened" Vector object
	 */
	public Vector toVector() {
		Vector vec = new Vector(rowCount * columnCount);
		double[] vecval = vec.getVector();
		int ii = 0;

		for (int i = 0; i < rowCount; i++)
			for (int j = 0; j < columnCount; j++)
				vecval[ii++] = this.data[i][j];

		return vec;
	}

	/**
	 * overridden toString()
	 */
	@Override
	public String toString() {
		StringBuilder sb = new StringBuilder();

		DecimalFormat df = new DecimalFormat("+#,##0.00;-#");
		df.setMinimumFractionDigits(8);

		for (int i = 0; i < rowCount; i++) {
			for (int j = 0; j < columnCount; j++) {
				sb.append(df.format(data[i][j]) + " ");
			}
			sb.append('\n');
		}

		return sb.toString();
	}

	/**
	 * clone this vector ;) Cloning is getting interesting! :)
	 * 
	 * @return the clone
	 */
	@Override
	public Object clone() {
		Matrix theCopy = new Matrix(rowCount, columnCount);

		for (int i = 0; i < rowCount; i++) {
			for (int j = 0; j < columnCount; j++) {
				theCopy.data[i][j] = data[i][j];
			}
		}

		theCopy.rowCount = rowCount;
		theCopy.columnCount = columnCount;

		return theCopy;
	}

	/**
	 * get a row vector
	 * 
	 * @param rowIndex
	 *            - the row index to obtain
	 * @return instance of Vector representing the requested row
	 */
	public Vector getRowVector(int rowIndex) {
		Vector row = new Vector(columnCount);

		double[] r = row.getVector();

		for (int i = 0; i < columnCount; i++)
			r[i] = data[rowIndex][i];

		return row;
	}

	/**
	 * get a column vector
	 * 
	 * @param colIndex
	 *            - the column index to obtain
	 * @return instance of Vector representing the requested column
	 */
	public Vector getColumnVector(int colIndex) {
		Vector column = new Vector(rowCount);

		double[] c = column.getVector();

		for (int i = 0; i < rowCount; i++)
			c[i] = data[i][colIndex];

		return column;
	}

	/**
	 * Getter for property rowCount.
	 * 
	 * @return Value of property rowCount.
	 */
	public int getRowDimension() {
		return this.rowCount;
	}

	/**
	 * Setter for property rowCount.
	 * 
	 * @param rowCount
	 *            New value of property rowCount.
	 */
	public void setRowCount(int rowCount) {
		this.rowCount = rowCount;
	}

	/**
	 * Getter for property columnCount.
	 * 
	 * @return Value of property columnCount.
	 */
	public int getColumnDimension() {
		return this.columnCount;
	}

	/**
	 * Setter for property columnCount.
	 * 
	 * @param columnCount
	 *            New value of property columnCount.
	 */
	public void setColumnCount(int columnCount) {
		this.columnCount = columnCount;
	}

	/**
	 * Is this upper triangular matrix?
	 * 
	 * @return true is upper triangular, false otherwise
	 */
	public boolean isUpperTriangular() {
		if (this.columnCount != this.rowCount)
			return false;

		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < i; j++) {
				if (data[i][j] != 0) {
					return false;
				}
			}
		}

		return true;
	}

	/**
	 * isSingular() - This method finds whether a matrix is singular or not i.e if
	 * pth row contains all zeros.
	 * 
	 * @param p
	 *            - The row to be checked for zero.
	 * @return boolean - true if the matrix is singular - false otherwise
	 */
	public boolean isSingular(int p) {
		for (int i = p; i <= data[p].length; i++) {
			if (data[p][i] != 0) {
				return false;
			}
		}

		return true;
	}

	/**
	 * isSingular() - This method finds whether a matrix is singular or not i.e if
	 * pth row contains all zeros.
	 * 
	 * @param p
	 *            - The row to be checked for zero.
	 * @param row
	 *            - modified row indices
	 * @return boolean - true if the matrix is singular - false otherwise
	 */
	public boolean isSingular(int p, int[] row) {
		if (row.length != data[row[p]].length)
			throw new UnsupportedOperationException("Size of row index and matrix do no match!");

		for (int i = p; i <= data[0].length - 1; i++) {
			if (data[row[p]][i] != 0) {
				return false;
			}
		}

		return true;
	}

	/**
	 * Returns the max norm of this matrix
	 * 
	 * @return max(abs(matrix))
	 */
	public double maxNorm() {
		double mx = Math.abs(data[0][0]);

		for (int i = 0; i < rowCount; i++)
			for (int j = 0; j < columnCount; j++)
				if (mx < Math.abs(data[i][j]))
					mx = Math.abs(data[i][j]);

		return mx;
	}

	/**
	 * Return the scaled instance of this matrix, each element multiplied by the
	 * scale factor
	 * 
	 * @param scaleFactor
	 *            the scale factor with which each element is multiplied
	 * @return the scaled instance of the current matrix
	 */
	public Matrix multiply(double scaleFactor) {
		Matrix res = new Matrix(rowCount, columnCount);

		for (int i = 0; i < rowCount; i++)
			for (int j = 0; j < columnCount; j++)
				res.data[i][j] = scaleFactor * data[i][j];

		return res;
	}
}