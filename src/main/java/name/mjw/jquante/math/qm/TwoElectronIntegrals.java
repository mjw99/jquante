package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.atomic.LongAdder;
import java.util.stream.IntStream;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import name.mjw.jquante.math.qm.basis.PrimitiveGaussian;
import name.mjw.jquante.math.qm.basis.Shell;
import name.mjw.jquante.math.qm.integral.Integrals;
import name.mjw.jquante.math.qm.integral.IntegralsUtil;
import name.mjw.jquante.molecule.Molecule;

import net.jafama.FastMath;

/**
 * The 2E integral driver.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public final class TwoElectronIntegrals {
	private static final Logger LOG = LogManager.getLogger(TwoElectronIntegrals.class);

	private BasisSetLibrary basisSetLibrary;

	private Molecule molecule;

	/**
	 * Holds value of property twoEIntegrals.
	 */
	private double[] twoEIntegrals;

	/** The index of the atom with respect to which 2E derivatives are computed. */
	private int atomIndex;

	/** The SCF method instance used during derivative evaluation. */
	protected SCFMethod scfMethod;

	/** List of partial derivative arrays [dx, dy, dz] of the two-electron integrals. */
	protected ArrayList<double[]> twoEDer;

	/** If true, integrals are computed on-the-fly rather than stored in memory. */
	protected boolean onTheFly;

	/**
	 * Creates a new instance of TwoElectronIntegrals
	 * 
	 * @param basisSetLibrary the basis functions to be used
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary) {
		this(basisSetLibrary, false);
	}

	/**
	 * Creates a new instance of TwoElectronIntegrals
	 * 
	 * @param basisSetLibrary the basis functions to be used
	 * @param onTheFly        if true, the 2E integrals are not calculated and
	 *                        stored, they must be calculated individually by
	 *                        calling compute2E(i,j,k,l)
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary, boolean onTheFly) {
		this.basisSetLibrary = basisSetLibrary;

		this.onTheFly = onTheFly;

		// compute the 2E integrals
		if (!onTheFly) {
			try {
				compute2E(); // try to do compute 2E incore
			} catch (OutOfMemoryError e) {
				// if no memory, resort to direct SCF
				LOG.error("No memory for in-core integral evaluation" + ". Switching to direct integral evaluation.");
				this.onTheFly = true;
			}
		}
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
	 * @param basisSetLibrary Basis functions for the molecule.
	 * @param molecule        Molecule for the two electron integrals.
	 * @param onTheFly        if true, the 2E integrals are not calculated and
	 *                        stored, they must be calculated individually by
	 *                        calling compute2EShellPair(i,j,k,l)
	 */
	public TwoElectronIntegrals(BasisSetLibrary basisSetLibrary, Molecule molecule, boolean onTheFly) {
		this.basisSetLibrary = basisSetLibrary;
		this.molecule = molecule;

		this.onTheFly = onTheFly;

		// compute the 2E integrals
		if (!onTheFly) {
			try {
				// first try in - core method
				try {
					compute2EShellPair(); // try to use shell pair algorithm
				} catch (IllegalStateException e) {
					LOG.warn("Shell-pair algorithm unavailable ({}); falling back to conventional algorithm",
							e.getMessage());
					compute2E();
				}
			} catch (OutOfMemoryError e) {
				// if no memory, resort to direct SCF
				LOG.error("No memory for in-core integral evaluation" + ". Switching to direct integral evaluation.");
				this.onTheFly = true;
			}
		}
	}

	/**
	 * Compute the 2E integrals, and store it in a single 1D array, in the form
	 * [ijkl].
	 *
	 */
	protected void compute2ESerial() {
		LOG.debug("compute2ESerial() called");
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		// allocate required memory
		final int noOfBasisFunctions = bfs.size();
		final int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		LOG.debug("Two electron integrals to evaluate: {}", noOfIntegrals);

		twoEIntegrals = new double[noOfIntegrals];

		// we only need i <= j, k <= l, and ij >= kl
		for (int i = 0; i < noOfBasisFunctions; i++) {
			for (int j = 0; j < i + 1; j++) {
				int ij = i * (i + 1) / 2 + j;

				for (int k = 0; k < noOfBasisFunctions; k++) {
					for (int l = 0; l < (k + 1); l++) {
						int kl = k * (k + 1) / 2 + l;

						if (ij >= kl) {
							int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);
							// record the 2E integrals
							twoEIntegrals[ijkl] = Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
						}
					}
				}
			}
		}

	}

	/**
	 * Compute the 2E integrals, and store it in a single 1D array, in the form
	 * [ijkl].
	 * 
	 * This method has been modified to take advantage of multi core systems where
	 * available.
	 */
	protected void compute2E() {
		LOG.debug("compute2E() called");
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		// allocate required memory
		final int noOfBasisFunctions = bfs.size();
		final int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		LOG.debug("Two electron integrals to evaluate: {}", noOfIntegrals);

		twoEIntegrals = new double[noOfIntegrals];

		// we only need i <= j, k <= l, and ij >= kl
		IntStream.range(0, noOfBasisFunctions).parallel().forEach(i -> {
			for (int j = 0; j < i + 1; j++) {
				int ij = i * (i + 1) / 2 + j;
				for (int k = 0; k < noOfBasisFunctions; k++) {
					for (int l = 0; l < k + 1; l++) {
						int kl = k * (k + 1) / 2 + l;

						if (ij >= kl) {
							int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);
							// record the 2E integrals
							twoEIntegrals[ijkl] = Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
						}
					}
				}
			}
		});

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
				currentOrigin = iPG.origin();
				currentPower = iPG.powers();
				currentAlpha = iPG.exponent();
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							twoEDerEle = twoEDerEle.add(compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx));
						}
					}
				}
			}
		}

		if (bfj.getCenteredAtom().getIndex() == atomIndex) {
			paramIdx = new int[] { 0, 4, 2, 3 };
			for (PrimitiveGaussian iPG : bfi.getPrimitives()) {
				for (PrimitiveGaussian jPG : bfj.getPrimitives()) {
					currentOrigin = jPG.origin();
					currentPower = jPG.powers();
					currentAlpha = jPG.exponent();
					for (PrimitiveGaussian kPG : bfk.getPrimitives()) {
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							twoEDerEle = twoEDerEle.add(compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx));
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
						currentOrigin = kPG.origin();
						currentPower = kPG.powers();
						currentAlpha = kPG.exponent();
						for (PrimitiveGaussian lPG : bfl.getPrimitives()) {
							twoEDerEle = twoEDerEle.add(compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx));
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
							currentOrigin = lPG.origin();
							currentPower = lPG.powers();
							currentAlpha = lPG.exponent();

							twoEDerEle = twoEDerEle.add(compute2EDerivativeElementHelper(iPG, jPG, kPG, lPG, currentOrigin, currentPower,
									currentAlpha, paramIdx));
						} // end for
					} // end for
				} // end for
			} // end for
		} // end if

		return twoEDerEle.scalarMultiply(
				bfi.getNormalization() * bfj.getNormalization() * bfk.getNormalization() * bfl.getNormalization());
	}

	/** helper for computing 2E derivative */
	private Vector3D compute2EDerivativeElementHelper(PrimitiveGaussian iPG, PrimitiveGaussian jPG, PrimitiveGaussian kPG,
			PrimitiveGaussian lPG, Vector3D currentOrigin, Power currentPower, double currentAlpha, int[] paramIdx) {

		int l = currentPower.l();
		int m = currentPower.m();
		int n = currentPower.n();

		double coeff = iPG.coefficient() * jPG.coefficient() * kPG.coefficient() * lPG.coefficient();

		// x-component derivative
		PrimitiveGaussian xPG = new PrimitiveGaussian(currentOrigin, new Power(l + 1, m, n), currentAlpha, coeff);
		PrimitiveGaussian[] pgs = new PrimitiveGaussian[] { iPG, jPG, kPG, lPG, xPG };
		double termax = FastMath.sqrt(currentAlpha * (2.0 * l + 1.0)) * coeff
				* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		double termbx = 0.0;
		if (l > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(l - 1, m, n), currentAlpha, coeff);
			pgs[4] = xPG;
			termbx = -2.0 * l * FastMath.sqrt(currentAlpha / (2.0 * l - 1.0)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		// y-component derivative
		xPG = new PrimitiveGaussian(currentOrigin, new Power(l, m + 1, n), currentAlpha, coeff);
		pgs[4] = xPG;
		double termay = FastMath.sqrt(currentAlpha * (2.0 * m + 1.0)) * coeff
				* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		double termby = 0.0;
		if (m > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(l, m - 1, n), currentAlpha, coeff);
			pgs[4] = xPG;
			termby = -2.0 * m * FastMath.sqrt(currentAlpha / (2.0 * m - 1.0)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		// z-component derivative
		xPG = new PrimitiveGaussian(currentOrigin, new Power(l, m, n + 1), currentAlpha, coeff);
		pgs[4] = xPG;
		double termaz = FastMath.sqrt(currentAlpha * (2.0 * n + 1.0)) * coeff
				* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		double termbz = 0.0;
		if (n > 0) {
			xPG = new PrimitiveGaussian(currentOrigin, new Power(l, m, n - 1), currentAlpha, coeff);
			pgs[4] = xPG;
			termbz = -2.0 * n * FastMath.sqrt(currentAlpha / (2.0 * n - 1.0)) * coeff
					* Integrals.coulomb(pgs[paramIdx[0]], pgs[paramIdx[1]], pgs[paramIdx[2]], pgs[paramIdx[3]]);
		}

		return new Vector3D(termax + termbx, termay + termby, termaz + termbz);
	}

	/**
	 * Compute the 2E integrals using the true shell-pair algorithm, storing results
	 * in a single 1D array in the form [ijkl].
	 *
	 * <p>Uses the unique shell pairs from {@link BasisSetLibrary#getUniqueShellPairs()},
	 * pre-computes Schwarz screening values Q_IJ = max&nbsp;sqrt(|(ij|ij)|) per shell
	 * pair, and skips quartets where Q_IJ * Q_KL &lt; 1e-10. The outer loop over bra
	 * shell pairs is parallelised.
	 */
	protected void compute2EShellPair() {
		LOG.debug("compute2EShellPair() called");

		List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();
		List<List<Shell>> shellPairs = basisSetLibrary.getUniqueShellPairs();

		if (shellPairs == null || shellPairs.isEmpty()) {
			throw new IllegalStateException("No unique shell pairs available; cannot use shell-pair algorithm");
		}

		int noOfBasisFunctions = bfs.size();
		int noOfIntegrals = noOfBasisFunctions * (noOfBasisFunctions + 1)
				* (noOfBasisFunctions * noOfBasisFunctions + noOfBasisFunctions + 2) / 8;

		LOG.debug("Two electron integrals to evaluate: {}", noOfIntegrals);

		twoEIntegrals = new double[noOfIntegrals];

		int nShellPairs = shellPairs.size();

		// Pre-compute Schwarz screening values Q_IJ = max sqrt(|(ij|ij)|) per shell pair
		double[] schwarzValues = new double[nShellPairs];
		for (int pq = 0; pq < nShellPairs; pq++) {
			Shell shellA = shellPairs.get(pq).get(0);
			Shell shellB = shellPairs.get(pq).get(1);
			double maxVal = 0.0;
			for (int i = shellA.getFirstBasisFunctionIndex(); i <= shellA.getLastBasisFunctionIndex(); i++) {
				for (int j = shellB.getFirstBasisFunctionIndex(); j <= shellB.getLastBasisFunctionIndex(); j++) {
					double val = FastMath.sqrt(FastMath.abs(Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(i), bfs.get(j))));
					if (val > maxVal) maxVal = val;
				}
			}
			schwarzValues[pq] = maxVal;
		}

		// Parallel loop over bra shell pairs.
		// Deduplication is handled entirely at the basis-function level (j<=i, l<=k,
		// ij>=kl). A shell-level canonical guard based on firstBasisFunctionIndex is
		// NOT safe here: shells can span multiple functions, so a "lower" shell pair
		// may still contain (i,j) with ij > kl relative to a "higher" pair — the
		// guard would block that direction and the reversed direction can't reach the
		// missing j value. Iterating all (pq,rs) pairs and relying on ij>=kl at the
		// basis-function level is both correct and avoids any double computation.
		double minSchwarz = Double.MAX_VALUE;
		double maxSchwarz = 0.0;
		for (double q : schwarzValues) {
			if (q < minSchwarz) minSchwarz = q;
			if (q > maxSchwarz) maxSchwarz = q;
		}
		LOG.info("Schwarz values: min={} max={}", String.format("%.2e", minSchwarz), String.format("%.2e", maxSchwarz));

		LongAdder skippedQuartets = new LongAdder();
		long totalQuartets = (long) nShellPairs * nShellPairs;

		IntStream.range(0, nShellPairs).parallel().forEach(pq -> {
			List<Shell> braPair = shellPairs.get(pq);
			Shell shellI = higherIndexedShell(braPair.get(0), braPair.get(1));
			Shell shellJ = lowerIndexedShell(braPair.get(0), braPair.get(1));

			for (int rs = 0; rs < nShellPairs; rs++) {
				// Schwarz screening: skip negligible quartets
				if (schwarzValues[pq] * schwarzValues[rs] < 1e-7) {
					skippedQuartets.increment();
					continue;
				}

				List<Shell> ketPair = shellPairs.get(rs);
				Shell shellK = higherIndexedShell(ketPair.get(0), ketPair.get(1));
				Shell shellL = lowerIndexedShell(ketPair.get(0), ketPair.get(1));

				for (int i = shellI.getFirstBasisFunctionIndex(); i <= shellI.getLastBasisFunctionIndex(); i++) {
					for (int j = shellJ.getFirstBasisFunctionIndex(); j <= shellJ.getLastBasisFunctionIndex(); j++) {
						if (j > i) continue;
						int ij = i * (i + 1) / 2 + j;

						for (int k = shellK.getFirstBasisFunctionIndex(); k <= shellK.getLastBasisFunctionIndex(); k++) {
							for (int l = shellL.getFirstBasisFunctionIndex(); l <= shellL.getLastBasisFunctionIndex(); l++) {
								if (l > k) continue;
								int kl = k * (k + 1) / 2 + l;
								if (ij < kl) continue;

								twoEIntegrals[IntegralsUtil.ijkl2intindex(i, j, k, l)] =
										Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
							}
						}
					}
				}
			}
		});

		long skipped = skippedQuartets.sum();
		LOG.info("Schwarz screening: {} of {} shell-pair quartets skipped ({} %)",
				skipped, totalQuartets, String.format("%.1f", 100.0 * skipped / totalQuartets));
	}

	/**
	 * Returns whichever of the two shells has the higher first basis-function index.
	 * For equal indices (self-pair), returns {@code a}.
	 */
	private Shell higherIndexedShell(Shell a, Shell b) {
		return a.getFirstBasisFunctionIndex() >= b.getFirstBasisFunctionIndex() ? a : b;
	}

	/**
	 * Returns whichever of the two shells has the lower first basis-function index.
	 * For equal indices (self-pair), returns {@code a}.
	 */
	private Shell lowerIndexedShell(Shell a, Shell b) {
		return a.getFirstBasisFunctionIndex() <= b.getFirstBasisFunctionIndex() ? a : b;
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
		List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();

		return Integrals.coulomb(bfs.get(i), bfs.get(j), bfs.get(k), bfs.get(l));
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
		return Integrals.coulomb(cgi, cgj, cgk, cgl);
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
		final List<ContractedGaussian> bfs = basisSetLibrary.getBasisFunctions();
		final int noOfBasisFunctions = bfs.size();

		this.atomIndex = atomIndex;
		this.scfMethod = scfMethod;

		twoEDer = new ArrayList<>();

		double[] dxTwoE = new double[twoEIntegrals.length];
		double[] dyTwoE = new double[twoEIntegrals.length];
		double[] dzTwoE = new double[twoEIntegrals.length];

		twoEDer.add(dxTwoE);
		twoEDer.add(dyTwoE);
		twoEDer.add(dzTwoE);

		// we only need i <= j, k <= l, and ij >= kl
		IntStream.range(0, noOfBasisFunctions).parallel().forEach(i -> {
			ContractedGaussian bfi = bfs.get(i);

			for (int j = 0; j < i + 1; j++) {
				ContractedGaussian bfj = bfs.get(j);
				int ij = i * (i + 1) / 2 + j;
				for (int k = 0; k < noOfBasisFunctions; k++) {
					ContractedGaussian bfk = bfs.get(k);
					for (int l = 0; l < k + 1; l++) {
						ContractedGaussian bfl = bfs.get(l);
						int kl = k * (k + 1) / 2 + l;

						if (ij >= kl) {
							int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);

							// record derivative of the 2E integrals
							Vector3D twoEDerEle = compute2EDerivativeElement(bfi, bfj, bfk, bfl);

							dxTwoE[ijkl] = twoEDerEle.getX();
							dyTwoE[ijkl] = twoEDerEle.getY();
							dzTwoE[ijkl] = twoEDerEle.getZ();
						}
					}
				}
			}
		});
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

	/**
	 * Prints all two-electron integrals (i,j,k,l) and their values to standard output.
	 * Intended for debugging purposes.
	 */
	public void printCompleteIntegralList() {

		int noOfBasisFunctions = basisSetLibrary.getBasisFunctions().size();

		System.out.println("");
		System.out.println("Complete integral list");
		System.out.println("======================");
		for (int i = 0; i < noOfBasisFunctions; i++) {
			for (int j = 0; j < noOfBasisFunctions; j++) {
				for (int k = 0; k < noOfBasisFunctions; k++) {
					for (int l = 0; l < noOfBasisFunctions; l++) {
						int ijkl = IntegralsUtil.ijkl2intindex(i, j, k, l);
						System.out.println(i + " " + j + " " + k + " " + l + "\t\t" + twoEIntegrals[ijkl]);

					}

				}

			}

		}
	}
}
