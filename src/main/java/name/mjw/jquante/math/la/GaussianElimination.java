package name.mjw.jquante.math.la;

import name.mjw.jquante.math.Matrix;
import name.mjw.jquante.math.la.exception.SingularMatrixException;

/**
 * The Gaussian elimination solver for A x = B
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class GaussianElimination extends LinearEquationSolver {

	private int[] row;
	private double[][] a;
	private int n;
	private int n1;

	@Override
	public void findSolution() throws SingularMatrixException {
		int N = matrixA.getRowCount();
		double[] x;
		Matrix m = new Matrix(N, N + 1);
		a = m.getMatrix();
		x = vectorX.getVector();

		double[][] ta = matrixA.getMatrix();
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				a[i][j] = ta[i][j];

		for (int i = 0; i < N; i++)
			a[i][N] = x[i];

		n = a.length - 1;
		n1 = a[0].length - 1;

		row = new int[a.length];

		// Initialise row vector
		for (int j = 1; j <= n; j++)
			row[j] = j;

		// first make the the matrix upper triangular
		upperTriangularize(false);

		// now use back substitution to get the solution
		// first solution from back :)
		x[n] = a[row[n]][n1] / a[row[n]][n];

		double sum;
		for (int k = (n - 1); k >= 0; k--) {
			sum = 0.0;
			for (int c = (k + 1); c <= n; c++) {
				sum += (a[row[k]][c] * x[c]); // compute next solutions
			}

			// the actual next solution
			x[k] = (a[row[k]][n1] - sum) / a[row[k]][k];

		}
	}

	/** Upper triangularize the matrix */
	private void upperTriangularize(boolean doOneScale) throws SingularMatrixException {
		// check if matrix is already upper triangular
		// this check may be removed in future
		if (matrixA.isUpperTriangular()) {
			return;
		}

		for (int p = 0; p <= n1; p++) {
			simplePivot(p); // apply simple pivoting
			oneScale(p, doOneScale); // apply simple 1-scaling
		}

		if (matrixA.isSingular(n, row)) {
			throw new SingularMatrixException(); // if singular
		}
	}

	/**
	 * This method does simple pivoting if the a[j][j]th element is zero. It first
	 * finds the a[i][j]th element which is numerically greater than a[j][j] and
	 * i&gt;j and then interchanges the ith and jth rows.
	 * 
	 * @param p
	 *            The pth iteration in Gaussian elemination.
	 * @throws SingularMatrixException
	 *             The matrix is singular.
	 */
	public void simplePivot(int p) throws SingularMatrixException {
		int temp = 0;

		if (p >= row.length) {
			return;
		}

		// if concerned row's first element is unity
		// or more do not pivot
		if (Math.abs(a[row[p]][p]) >= 1)
			return;

		for (int k = (p + 1); k <= n; k++) {
			// check for keeping things near unity
			if ((Math.abs(a[row[k]][p]) - 1) < (Math.abs(a[row[p]][p]) - 1)) {
				// switch the indices to represent a
				// row interchange so as to avoid the overhead
				// of actually moving the row elements :)
				temp = row[p];
				row[p] = row[k];
				row[k] = temp;
			}
		}

		// check if singular
		if (matrixA.isSingular(p, row)) {
			throw new SingularMatrixException(); // if singular
		}
	}

	/**
	 * This method makes the a[j][j]th entry at the jth iteration close to unity by
	 * dividing each element of the jth row by a[j][j].
	 * 
	 * @param p
	 *            - The jth iteration in Gaussian elimation. boolean scale - Scale
	 *            to one or not
	 */
	private void oneScale(int p, boolean scale) {
		double m;

		if (p >= row.length) {
			return;
		}

		for (int k = (p + 1); k <= n; k++) {
			// compute the multiplicity factor
			m = a[row[k]][p] / a[row[p]][p];

			for (int c = (p + 1); c <= n1; c++) {
				a[row[k]][c] -= (m * a[row[p]][c]);
			}

			a[row[k]][p] = m;
		}

		if (scale) {
			a[row[p]][p] = 1.0; // diagonal element is now 1
		}

		for (int c = (p + 1); c <= n; c++) {
			a[row[c]][p] = 0.0; // all entries below diagonal are 0
		}
	}
}
