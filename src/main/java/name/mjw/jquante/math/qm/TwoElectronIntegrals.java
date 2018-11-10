package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.qm.basis.CompactContractedGaussian;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import name.mjw.jquante.math.qm.basis.PrimitiveGaussian;
import name.mjw.jquante.math.qm.integral.Integrals;
import name.mjw.jquante.math.qm.integral.IntegralsUtil;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.parallel.AbstractSimpleParallelTask;
import name.mjw.jquante.parallel.SimpleParallelTask;
import name.mjw.jquante.parallel.SimpleParallelTaskExecuter;
import net.jafama.FastMath;

/**
 * The 2E integral driver.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class TwoElectronIntegrals {
	private static final Logger LOG = LogManager.getLogger(TwoElectronIntegrals.class);

	private BasisFunctions basisFunctions;

	private Molecule molecule;

	/**
	 * Holds value of property twoEIntegrals.
	 */
	private double[] twoEIntegrals;

	private int atomIndex;
	protected SCFMethod scfMethod;

	protected ArrayList<double[]> twoEDer;

	/** Calculate integrals on the fly instead of storing them in an array */
	protected boolean onTheFly;

	/**
	 * Creates a new instance of TwoElectronIntegrals
	 * 
	 * @param basisFunctions the basis functions to be used
	 */
	public TwoElectronIntegrals(BasisFunctions basisFunctions) {
		this(basisFunctions, false);
	}

	/**
	 * Creates a new instance of TwoElectronIntegrals
	 * 
	 * @param basisFunctions the basis functions to be used
	 * @param onTheFly       if true, the 2E integrals are not calculated and
	 *                       stored, they must be calculated individually by calling
	 *                       compute2E(i,j,k,l)
	 */
	public TwoElectronIntegrals(BasisFunctions basisFunctions, boolean onTheFly) {
		this.basisFunctions = basisFunctions;

		this.onTheFly = onTheFly;

		// compute the 2E integrals
		if (!onTheFly) {
			try {
				compute2E(); // try to do compute 2E incore
			} catch (OutOfMemoryError e) {
				// if no memory, resort to direct SCF
				LOG.error("No memory for in-core integral evaluation" + ". Switching to direct integral evaluation.");
				this.onTheFly = true;
			} // end of try catch block
		} // end if
	}

	/**
	 * Creates a new instance of TwoElectronIntegrals, and computes two electron
	 * integrals using shell-pair method rather than the conventional loop over all
	 * basis functions approach. Note that, this requires the Molecule object to be
	 * passed as an argument for it to function correctly. Note that the Atom
	 * objects which are a part of Molecule object must be specially configured to
	 * have a user defined property called "basisFunctions" that is essentially a
	 * list of ContractedGaussians that are centered on the atom. If this breaks,
	 * then the code will silently default to using the conventional way of
	 * computing the two electron integrals.
	 * 
	 * @param basisFunctions Basis functions for the molecule.
	 * @param molecule       Molecule for the two electron integrals.
	 * @param onTheFly       if true, the 2E integrals are not calculated and
	 *                       stored, they must be calculated individually by calling
	 *                       compute2EShellPair(i,j,k,l)
	 */
	public TwoElectronIntegrals(BasisFunctions basisFunctions, Molecule molecule, boolean onTheFly) {
		this.basisFunctions = basisFunctions;
		this.molecule = molecule;

		this.onTheFly = onTheFly;

		// compute the 2E integrals
		if (!onTheFly) {
			try {
				// first try in - core method
				try {
					compute2EShellPair(); // try to use shell pair algorithm
				} catch (Exception ignored) {
					LOG.error("WARNING: Using a slower algorithm" + " to compute two electron integrals. Error is: "
							+ ignored.toString());
					LOG.error(ignored);

					System.out.println("Calling compute2E()");
					compute2E(); // if it doesn't succeed then fall back to
									// normal
				} // end if
			} catch (OutOfMemoryError e) {
                                 System.out.println("Switching to direct integral evaluation");
				// if no memory, resort to direct SCF
				LOG.error("No memory for in-core integral evaluation" + ". Switching to direct integral evaluation.");
				this.onTheFly = true;
			}
		}
	}

	/**
	 * compute the 2E integrals, and store it in a single 1D array, in the form
	 * [ijkl].
	 * 
	 * This method has been modified to take advantage of multi core systems where
	 * available.
	 */
	protected void compute2E() {
		final List<ContractedGaussian> bfs = basisFunctions.getBasisFunctions();

		// allocate required memory
		final int noOfBasisFunctions = bfs.size();
		final int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		twoEIntegrals = new double[noOfIntegrals];

		final SimpleParallelTaskExecuter pTaskExecuter = new SimpleParallelTaskExecuter();

		final TwoElectronIntegralEvaluaterThread tThread = new TwoElectronIntegralEvaluaterThread();
		tThread.setTaskName("TwoElectronIntegralEvaluater Thread");
		tThread.setTotalItems(noOfBasisFunctions);

		pTaskExecuter.execute(tThread);
	}

	/** Compute a single 2E derivative element */
	private Vector3D compute2EDerivativeElement(ContractedGaussian bfi, ContractedGaussian bfj, ContractedGaussian bfk,
			ContractedGaussian bfl) {

		Vector3D twoEDerEle = new Vector3D(0, 0, 0);
		int[] paramIdx;
		Power currentPower;
		Vector3D currentOrigin;
		double currentAlpha;

		if (bfi.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 4, 1, 2, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				currentOrigin = iPG.getOrigin();
				currentPower = iPG.getPowers();
				currentAlpha = iPG.getExponent();
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						}
					}
				}
			}
		}

		if (bfj.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 4, 2, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					currentOrigin = jPG.getOrigin();
					currentPower = jPG.getPowers();
					currentAlpha = jPG.getExponent();
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		if (bfk.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 1, 4, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						currentOrigin = kPG.getOrigin();
						currentPower = kPG.getPowers();
						currentAlpha = kPG.getExponent();
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		if (bfl.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 1, 2, 4 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							currentOrigin = lPG.getOrigin();
							currentPower = lPG.getPowers();
							currentAlpha = lPG.getExponent();

							compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx, twoEDerEle);
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		return twoEDerEle;
	}

	/** helper for computing 2E derivative */
	private void compute2EDerivativeElementHelper(PrimitiveGaussian iPG, PrimitiveGaussian jPG, PrimitiveGaussian kPG,
			PrimitiveGaussian lPG, Vector3D currentOrigin, Power currentPower, double currentAlpha, int[] paramIdx,
			Vector3D derEle) {

		int l = currentPower.getL();
		int m = currentPower.getM();
		int n = currentPower.getN();

		double coeff = iPG.getCoefficient() * jPG.getCoefficient() * kPG.getCoefficient() * lPG.getCoefficient();

		PrimitiveGaussian xPG = new PrimitiveGaussian(currentOrigin, new Power(l + 1, m, n), currentAlpha, coeff);
		xPG.normalize();

		PrimitiveGaussian[] pgs = new PrimitiveGaussian[] { iPG, jPG, kPG, lPG, xPG };

		double terma = FastMath.sqrt(currentAlpha * (2.0 * l + 1.0)) * coeff
				* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		double termbx = 0.0;
		double termby = 0.0;
		double termbz = 0.0;

		if (l > 0) {
			xPG.setPowers(new Power(l - 1, m, n));
			xPG.normalize();
			termbx = -2.0 * l * FastMath.sqrt(currentAlpha / (2. * l - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		} // end if

		if (m > 0) {
			xPG.setPowers(new Power(l, m - 1, n));
			xPG.normalize();
			termby = -2.0 * m * FastMath.sqrt(currentAlpha / (2. * m - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		} // end if

		if (n > 0) {
			xPG.setPowers(new Power(l, m, n - 1));
			xPG.normalize();
			termbz = -2.0 * n * FastMath.sqrt(currentAlpha / (2. * n - 1)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		} // end if

		derEle = new Vector3D(derEle.getX() + terma + termbx, derEle.getY() + terma + termby,
				derEle.getZ() + terma + termbz);
	}

	/**
	 * compute the 2E integrals using shell pair based method, and store it in a
	 * single 1D array, in the form [ijkl].
	 * 
	 * This method has been modified to take advantage of multi core systems where
	 * available.
	 */
	protected void compute2EShellPair() {
		// TODO : parallel

		// System.out.println("compute2EShellPair()");
		ArrayList<ContractedGaussian> bfs = basisFunctions.getBasisFunctions();

		// allocate required memory
		int noOfBasisFunctions = bfs.size();
		int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		twoEIntegrals = new double[noOfIntegrals];

		int noOfAtoms = molecule.getNumberOfAtoms();

		int naFunc;
		int nbFunc;
		int ncFunc;
		int ndFunc;
		int twoEIndx;

		ArrayList<ContractedGaussian> aFunc;
		ArrayList<ContractedGaussian> bFunc;
		ArrayList<ContractedGaussian> cFunc;
		ArrayList<ContractedGaussian> dFunc;
		ContractedGaussian iaFunc;
		ContractedGaussian jbFunc;
		ContractedGaussian kcFunc;
		ContractedGaussian ldFunc;

		// center a
		for (int a = 0; a < noOfAtoms; a++) {
			aFunc = (ArrayList<ContractedGaussian>) molecule.getAtom(a).getUserDefinedAtomProperty("basisFunctions")
					.getValue();
			naFunc = aFunc.size();

			// basis functions on a
			for (int i = 0; i < naFunc; i++) {
				iaFunc = aFunc.get(i);

				// center b
				for (int b = 0; b <= a; b++) {
					bFunc = (ArrayList<ContractedGaussian>) molecule.getAtom(b)
							.getUserDefinedAtomProperty("basisFunctions").getValue();
					nbFunc = (b < a) ? bFunc.size() : i + 1;
					// basis functions on b
					for (int j = 0; j < nbFunc; j++) {
						jbFunc = bFunc.get(j);

						// center c
						for (int c = 0; c < noOfAtoms; c++) {
							cFunc = (ArrayList<ContractedGaussian>) molecule.getAtom(c)
									.getUserDefinedAtomProperty("basisFunctions").getValue();
							ncFunc = cFunc.size();
							// basis functions on c
							for (int k = 0; k < ncFunc; k++) {
								kcFunc = cFunc.get(k);

								// center d
								for (int d = 0; d <= c; d++) {
									dFunc = (ArrayList<ContractedGaussian>) molecule.getAtom(d)
											.getUserDefinedAtomProperty("basisFunctions").getValue();
									ndFunc = (d < c) ? dFunc.size() : k + 1;
									// basis functions on d
									for (int l = 0; l < ndFunc; l++) {
										ldFunc = dFunc.get(l);

										twoEIndx = IntegralsUtil.ijkl2intindex(iaFunc.getIndex(), jbFunc.getIndex(),
												kcFunc.getIndex(), ldFunc.getIndex());

										twoEIntegrals[twoEIndx] = compute2E(iaFunc, jbFunc, kcFunc, ldFunc);

										//System.out.println(compute2E(iaFunc, jbFunc, kcFunc, ldFunc));
									} // end for (l)
								} // end for (d)
							} // end for (k)
						} // end for (c)
					} // end for (j)
				} // end for (b)
			} // end for (i)
		} // end for (a)
	}

	/**
	 * 
	 * Compute an integral centered at &lt;ij&#x7C;kl&gt;
	 * 
	 * @param i Index of contracted Gaussian function i.
	 * @param j Index of contracted Gaussian function j.
	 * @param k Index of contracted Gaussian function k.
	 * @param l Index of contracted Gaussian function l.
	 * @return the value of two integral electron
	 */
	public double compute2E(int i, int j, int k, int l) {
		// if we really need to compute, go ahead!
		ArrayList<ContractedGaussian> bfs = basisFunctions.getBasisFunctions();

		return Integrals.coulomb(bfs.get(i).toCompactContractedGaussian(), 
						bfs.get(j).toCompactContractedGaussian(), 
						bfs.get(k).toCompactContractedGaussian(), 
						bfs.get(l).toCompactContractedGaussian());
	}

	/**
	 * Compute an integral centered at &lt;ij&#x7C;kl&gt;
	 * 
	 * @param cgi Contracted Gaussian function i.
	 * @param cgj Contracted Gaussian function j.
	 * @param cgk Contracted Gaussian function k.
	 * @param cgl Contracted Gaussian function l.
	 * @return the value of two integral electron
	 */
	public double compute2E(ContractedGaussian cgi, ContractedGaussian cgj, ContractedGaussian cgk,
			ContractedGaussian cgl) {

		return Integrals.coulomb(cgi.toCompactContractedGaussian(), cgj.toCompactContractedGaussian(),
				cgk.toCompactContractedGaussian(), cgl.toCompactContractedGaussian());
	}

	/**
	 * Getter for property twoEIntegrals.
	 * 
	 * @return Value of property twoEIntegrals.
	 */
	public double[] getTwoEIntegrals() {
		return this.twoEIntegrals;
	}

	/**
	 * Setter for property twoEIntegrals.
	 * 
	 * @param twoEIntegrals New value of property twoEIntegrals.
	 */
	public void setTwoEIntegrals(double[] twoEIntegrals) {
		this.twoEIntegrals = twoEIntegrals;
	}

	/**
	 * Return the derivatives of the two electron integrals.
	 * 
	 * @param atomIndex the reference atom index
	 * @param scfMethod the reference SCF method
	 */
	public void compute2EDerivatives(int atomIndex, SCFMethod scfMethod) {
		this.atomIndex = atomIndex;
		this.scfMethod = scfMethod;

		twoEDer = new ArrayList<>();

		double[] dxTwoE = new double[twoEIntegrals.length];
		double[] dyTwoE = new double[twoEIntegrals.length];
		double[] dzTwoE = new double[twoEIntegrals.length];

		twoEDer.add(dxTwoE);
		twoEDer.add(dyTwoE);
		twoEDer.add(dzTwoE);

		SimpleParallelTaskExecuter pTaskExecuter = new SimpleParallelTaskExecuter();

		TwoElectronIntegralDerivativeEvaluaterThread tThread = new TwoElectronIntegralDerivativeEvaluaterThread();
		tThread.setTaskName("TwoElectronIntegralDerivativeEvaluater Thread");
		tThread.setTotalItems(basisFunctions.getBasisFunctions().size());

		pTaskExecuter.execute(tThread);
	}

	/**
	 * The currently computed list of two electron derivatives
	 * 
	 * @return partial derivatives of 2E integrals, computed in previous call to
	 *         compute2EDerivatives()
	 */
	public ArrayList<double[]> getTwoEDer() {
		return twoEDer;
	}

	/**
	 * Class encapsulating the way to compute 2E electrons in a way useful for
	 * utilizing multi core (processor) systems.
	 */
	protected class TwoElectronIntegralEvaluaterThread extends AbstractSimpleParallelTask {

		private int startBasisFunction;
		private int endBasisFunction;
		private ArrayList<ContractedGaussian> bfs;

		public TwoElectronIntegralEvaluaterThread() {
		}

		public TwoElectronIntegralEvaluaterThread(final int startBasisFunction, final int endBasisFunction,
				final ArrayList<ContractedGaussian> bfs) {
			this.startBasisFunction = startBasisFunction;
			this.endBasisFunction = endBasisFunction;

			this.bfs = bfs;

			setTaskName("TwoElectronIntegralEvaluater Thread");
		}

		/**
		 * Overridden run()
		 */
		@Override
		public void run() {
			compute2E(startBasisFunction, endBasisFunction, bfs);
		}

		/** Overridden init() */
		@Override
		public SimpleParallelTask init(final int startItem, final int endItem) {
			return new TwoElectronIntegralEvaluaterThread(startItem, endItem, basisFunctions.getBasisFunctions());
		}

		/**
		 * Actually compute the 2E integrals
		 */
		private void compute2E(final int startBasisFunction, final int endBasisFunction,
				final ArrayList<ContractedGaussian> bfs) {

			final int noOfBasisFunctions = bfs.size();

			int i;
			int j;
			int k;
			int l;

			int ij;
			int kl;

			int ijkl;

			final ArrayList<CompactContractedGaussian> cbfs = new ArrayList<>();

			for (ContractedGaussian item : bfs) {
				cbfs.add(item.toCompactContractedGaussian());
			}

			// we only need i <= j, k <= l, and ij >= kl
			for (i = startBasisFunction; i < endBasisFunction; i++) {
				for (j = 0; j < (i + 1); j++) {
					ij = i * (i + 1) / 2 + j;
					for (k = 0; k < noOfBasisFunctions; k++) {
						for (l = 0; l < (k + 1); l++) {
							kl = k * (k + 1) / 2 + l;
							if (ij >= kl) {
								ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);

								// record the 2E integrals
								twoEIntegrals[ijkl] = Integrals.coulomb(cbfs.get(i), cbfs.get(j), cbfs.get(k),
										cbfs.get(l));

							}
						}
					}
				}
			}

		}

	}

	/**
	 * Class encapsulating the way to compute 2E electrons in a way useful for
	 * utilizing multi core (processor) systems.
	 */
	protected class TwoElectronIntegralDerivativeEvaluaterThread extends AbstractSimpleParallelTask {

		private int startBasisFunction;
		private int endBasisFunction;
		private List<ContractedGaussian> bfs;

		public TwoElectronIntegralDerivativeEvaluaterThread() {
		}

		public TwoElectronIntegralDerivativeEvaluaterThread(int startBasisFunction, int endBasisFunction,
				List<ContractedGaussian> bfs) {
			this.startBasisFunction = startBasisFunction;
			this.endBasisFunction = endBasisFunction;

			this.bfs = bfs;

			setTaskName("TwoElectronIntegralDerivativeEvaluater Thread");
		}

		/**
		 * Actually compute the 2E integrals derivatives
		 */
		private void compute2EDerivative(int startBasisFunction, int endBasisFunction, List<ContractedGaussian> bfs) {

			int i;
			int j;
			int k;
			int l;
			int ij;
			int kl;
			int ijkl;
			int noOfBasisFunctions = bfs.size();

			ContractedGaussian bfi;
			ContractedGaussian bfj;
			ContractedGaussian bfk;
			ContractedGaussian bfl;

			double[] dxTwoE = twoEDer.get(0);
			double[] dyTwoE = twoEDer.get(1);
			double[] dzTwoE = twoEDer.get(2);

			// we only need i <= j, k <= l, and ij >= kl
			for (i = startBasisFunction; i < endBasisFunction; i++) {
				bfi = bfs.get(i);

				for (j = 0; j < (i + 1); j++) {
					bfj = bfs.get(j);
					ij = i * (i + 1) / 2 + j;

					for (k = 0; k < noOfBasisFunctions; k++) {
						bfk = bfs.get(k);

						for (l = 0; l < (k + 1); l++) {
							bfl = bfs.get(l);

							kl = k * (k + 1) / 2 + l;
							if (ij >= kl) {
								ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);

								// record derivative of the 2E integrals
								Vector3D twoEDerEle = compute2EDerivativeElement(bfi, bfj, bfk, bfl);

								dxTwoE[ijkl] = twoEDerEle.getX();
								dyTwoE[ijkl] = twoEDerEle.getY();
								dzTwoE[ijkl] = twoEDerEle.getZ();
							} // end if
						} // end l loop
					} // end k loop
				} // end of j loop
			} // end of i loop
		}

		/**
		 * Overridden run()
		 */
		@Override
		public void run() {
			compute2EDerivative(startBasisFunction, endBasisFunction, bfs);
		}

		/** Overridden init() */
		@Override
		public SimpleParallelTask init(int startItem, int endItem) {
			return new TwoElectronIntegralDerivativeEvaluaterThread(startItem, endItem,
					basisFunctions.getBasisFunctions());
		}
	} // end of class TwoElectronIntegralEvaluaterThread

	/**
	 * Get the value of onTheFly
	 * 
	 * @return the value of onTheFly
	 */
	public boolean isOnTheFly() {
		return onTheFly;
	}

	/**
	 * Set the value of onTheFly
	 * 
	 * @param onTheFly new value of onTheFly
	 */
	public void setOnTheFly(boolean onTheFly) {
		this.onTheFly = onTheFly;
	}

}
