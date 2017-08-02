package name.mjw.jquante.math.qm;

import name.mjw.jquante.math.Matrix;

/**
 * Represents the Density P matrix or the charge order bond density matrix (DM).
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class Density extends Matrix {

	/**
	 * Creates a new instance of NxM Matrix
	 * 
	 * @param n
	 *            the first dimension
	 * @param m
	 *            the second dimension
	 */
	public Density(int n, int m) {
		super(n, m);
	}

	/**
	 * Creates a new instance of square (NxN) Matrix
	 * 
	 * @param n
	 *            the dimension
	 */
	public Density(int n) {
		super(n, n);
	}

	/**
	 * Creates a new instance of Matrix, based on already allocated 2D array
	 * 
	 * @param a
	 *            the 2D array
	 */
	public Density(double[][] a) {
		super(a);
	}

	/**
	 * Compute the density matrix
	 * 
	 * @param scfMethod
	 *            The SCF method for which the density is to be formed
	 * @param guessInitialDM
	 *            Do a guess for DM?
	 * @param densityGuesser
	 *            Is there any special way to make the guess?
	 * @param noOfOccupancies
	 *            Occupied electrons
	 * @param mos
	 *            the MolecularOrbitals needed to actually compute the DM
	 */
	public void compute(SCFMethod scfMethod, boolean guessInitialDM,
			DensityGuesser densityGuesser, int noOfOccupancies,
			MolecularOrbitals mos) {
		if (guessInitialDM) {
			if (densityGuesser != null) {
				this.setMatrix(densityGuesser.guessDM(scfMethod).getMatrix());
				return;
			}
		}

		// else construct it from the MOs .. C*C'
		Matrix dVector = new Matrix(noOfOccupancies, mos.getRowCount());
		double[][] d = dVector.getMatrix();
		double[][] c = mos.getMatrix();

		int i, j;
		for (i = 0; i < noOfOccupancies; i++)
			for (j = 0; j < c.length; j++)
				d[i][j] = c[i][j];

		this.setMatrix(dVector.transpose().mul(dVector).getMatrix());
	}
}
