package name.mjw.jquante.math.qm;

import org.hipparchus.linear.ArrayRealVector;
import org.hipparchus.linear.EigenDecomposition;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.linear.RealVector;

/**
 * Reorders commons-math's EigenDecomposition eigensystem, sorting by the
 * smallest eigenvalue first.
 *
 * @author mjw
 *
 */
class SortedEigenDecomposition {

	double[] realEigenvalues;
	RealVector[] eigenvectors;

	final int n;

	public SortedEigenDecomposition(final EigenDecomposition eig) {
		n = eig.getRealEigenvalues().length;

		realEigenvalues = eig.getRealEigenvalues();
		eigenvectors = new ArrayRealVector[n];

		for (int k = 0; k < n; ++k) {
			eigenvectors[k] = eig.getEigenvector(k);
		}

		RealVector tmpvector;
		for (int i = 0; i < n; i++) {

			int k = i;
			double p = realEigenvalues[i];

			for (int j = i + 1; j < n; j++) {
				if (realEigenvalues[j] < p) {
					k = j;
					p = realEigenvalues[j];
				}
			}

			if (k != i) {
				realEigenvalues[k] = realEigenvalues[i];
				realEigenvalues[i] = p;

				tmpvector = eigenvectors[i];
				eigenvectors[i] = eigenvectors[k];
				eigenvectors[k] = tmpvector;
			}
		}

	}

	RealMatrix getVT() {
		RealMatrix vt = MatrixUtils.createRealMatrix(n, n);
		for (int k = 0; k < n; ++k) {
			vt.setRowVector(k, eigenvectors[k]);
		}

		return vt;
	}

	double[] getRealEigenvalues() {
		return realEigenvalues;
	}

}
