package name.mjw.jquante.math.qm;

import org.apache.log4j.Logger;

import name.mjw.jquante.math.Matrix;
import name.mjw.jquante.math.Vector3D;
import name.mjw.jquante.math.qm.event.SCFEvent;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.UserDefinedAtomProperty;

/**
 * Implements the Hartree-Fock (HF) SCF method for single point energy
 * evaluation of a molecule.
 *
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 *
 * @see https://en.wikipedia.org/wiki/Hartree%E2%80%93Fock_method
 * @see http://sirius.chem.vt.edu/~crawdad/programming/scf.pdf
 * @see http://sirius.chem.vt.edu/wiki/doku.php?id=crawdad:programming:project3
 */
public class HartreeFockSCFMethod extends SCFMethod {

	/**
	 * Logger object
	 */
	private static final Logger LOG = Logger
			.getLogger(HartreeFockSCFMethod.class);

	/**
	 * Represents an event in an SCF cycle
	 */
	protected SCFEvent scfEvent;

	/**
	 * SCF Type enumeration.
	 */
	private SCFType scfType;

	private boolean isDerivativeComputed = false;

	/**
	 * Creates a new instance of HartreeFockSCFMethod
	 * 
	 * @param molecule
	 *            The molecule under consideration.
	 * @param oneEI
	 *            The one electron integrals of the system
	 * @param twoEI
	 *            The two electron integrals of the system
	 */
	public HartreeFockSCFMethod(Molecule molecule, OneElectronIntegrals oneEI,
			TwoElectronIntegrals twoEI) {
		this(molecule, oneEI, twoEI, SCFType.HARTREE_FOCK);
	}

	/**
	 * Creates a new instance of HartreeFockSCFMethod
	 * 
	 * @param molecule
	 *            The molecule under consideration.
	 * @param oneEI
	 *            The one electron integrals of the system
	 * @param twoEI
	 *            The two electron integrals of the system
	 * @param scfType
	 *            Type of SCF Calculation.
	 */
	public HartreeFockSCFMethod(Molecule molecule, OneElectronIntegrals oneEI,
			TwoElectronIntegrals twoEI, SCFType scfType) {
		super(molecule, oneEI, twoEI);

		scfEvent = new SCFEvent(this);
		this.scfType = scfType;
	}

	/**
	 * Perform the SCF optimization of the molecular wave function until the
	 * energy converges.
	 */
	@Override
	public void scf() {
		// check first if closed shell run?
		int noOfElectrons = molecule.getNumberOfElectrons();
		int noOfOccupancies = noOfElectrons / 2;

		if (noOfElectrons % 2 != 0) {
			throw new UnsupportedOperationException("Open shell systems are"
					+ " not currently supported.");
		}

		Overlap overlap = oneEI.getOverlap();
		LOG.debug("Initial overlap\n" + overlap);

		HCore hCore = oneEI.getHCore();
		LOG.debug("Initial hCore\n" + hCore);

		boolean converged = false;
		double oldEnergy = 0.0;
		double nuclearEnergy = nuclearEnergy();
		double eOne;
		double eTwo;

		// init memory for the matrices
		gMatrix = new GMatrix(hCore.getRowCount());
		mos = new MolecularOrbitals(hCore.getRowCount());
		density = new Density(hCore.getRowCount());
		fock = new Fock(hCore.getRowCount());

		// compute initial MOs
		mos.compute(hCore, overlap);
		LOG.debug("Initial computed MO coefficient matrix as: \n" + mos);

		FockExtrapolator diis = new DIISFockExtrapolator();

		diis.init();

		LOG.debug("Initial density matrix \n" + density);

		// start the SCF cycle
		for (scfIteration = 0; scfIteration < maxIteration; scfIteration++) {

			LOG.debug("");
			LOG.debug("SCF iteration: " + scfIteration);
			// make or guess density
			density.compute(this, guessInitialDM && (scfIteration == 0),
					densityGuesser, noOfOccupancies, mos);

			// make the G matrix
			gMatrix.compute(scfType, twoEI, density);

			// make fock matrix
			fock.compute(hCore, gMatrix);

			// apply DIIS
			fock = diis.next(fock, overlap, density);

			// compute the new MOs
			mos.compute(fock, overlap);

			// compute the total energy at this point
			eOne = density.mul(hCore).trace();
			eTwo = density.mul(fock).trace();

			energy = eOne + eTwo + nuclearEnergy;

			LOG.debug("Energy is : " + energy + "\tdelta_E: "
					+ (energy - oldEnergy));

			// fire the SCF event notification
			scfEvent.setType(SCFEvent.INFO_EVENT);
			scfEvent.setCurrentIteration(scfIteration);
			scfEvent.setCurrentEnergy(energy);
			fireSCFEventListenerScfEventOccured(scfEvent);

			// check for convergence
			if (Math.abs(energy - oldEnergy) < energyTolerance) {
				converged = true;
				scfEvent.setType(SCFEvent.CONVERGED_EVENT);
				scfEvent.setCurrentIteration(scfIteration);
				scfEvent.setCurrentEnergy(energy);
				fireSCFEventListenerScfEventOccured(scfEvent);
				break;
			}

			oldEnergy = energy;
		} // end of SCF iteration

		// not converged? then inform so...
		if (!converged) {
			scfEvent.setType(SCFEvent.FAILED_CONVERGENCE_EVENT);
			scfEvent.setCurrentIteration(scfIteration);
			scfEvent.setCurrentEnergy(energy);
			fireSCFEventListenerScfEventOccured(scfEvent);
		}
	}

	/**
	 * This gradient (or Force) calculation is based on Appendix C of Mordern
	 * Quantum Chemistry by Szabo and Ostland, which describes computing
	 * analytic gradients and geometry optimization.
	 */
	private void computeForce() {
		Force hfForce = new HartreeFockForce();

		for (int i = 0; i < molecule.getNumberOfAtoms(); i++) {
			Atom atom = molecule.getAtom(i);
			Vector3D force = hfForce.computeForce(atom.getIndex(), this);

			UserDefinedAtomProperty atmForce = atom
					.getUserDefinedAtomProperty("force");
			if (atmForce == null) {
				atmForce = new UserDefinedAtomProperty("force", force);
				atom.addUserDefinedAtomProperty(atmForce);
			}

			atmForce.setValue(force);
		}

		isDerivativeComputed = true;
	}

	/**
	 * Evaluate the function with 'n' variables
	 * 
	 * @param variables
	 *            an array of variables
	 * @return a double value evaluating F(x1, x2, ...)
	 */
	@Override
	public double evaluate(double[] variables) {
		molecule.resetAtomCoordinates(variables, false);

		// perform scf, and return energy
		oneEI.compute1E();
		twoEI.compute2E();
		scf();
		isDerivativeComputed = false;
		return getEnergy();
	}

	/**
	 * For a method that encapsulates its own set of "base" variables (e.g. an
	 * initial set of atom positions), resets the base variables to new values.
	 * 
	 * @param variables
	 *            the new set of base variables
	 */
	@Override
	public void resetVariables(double[] variables) {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Getter for property derivativeAvailable.
	 * 
	 * @return Value of property derivativeAvailable.
	 */
	@Override
	public boolean isDerivativeAvailable() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Getter for property hessianAvailable.
	 * 
	 * @return Value of property hessianAvailable.
	 */
	@Override
	public boolean isHessianAvailable() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Getter for property derivatives.
	 * 
	 * @return Value of property derivatives.
	 */
	@Override
	public double[] getDerivatives() {
		if (!isDerivativeComputed) {
			computeForce();
		}

		// unpack and return the forces
		double[] forces = new double[molecule.getNumberOfAtoms() * 3];
		int ii = 0;
		for (int i = 0; i < molecule.getNumberOfAtoms(); i++) {
			UserDefinedAtomProperty atmForce = molecule.getAtom(i)
					.getUserDefinedAtomProperty("force");

			Vector3D force = (Vector3D) atmForce.getValue();
			forces[ii] = force.getI();
			forces[ii + 1] = force.getJ();
			forces[ii + 2] = force.getK();

			ii += 3;
		}

		return forces;
	}

	/**
	 * Return the max norm of the derivatives, if available
	 * 
	 * @return max norm of all the derivatives
	 */
	@Override
	public double getMaxNormOfDerivatives() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Return the RMS of the derivatives, if available
	 * 
	 * @return max RMS all the derivatives
	 */
	@Override
	public double getRMSOfDerivatives() {
		throw new UnsupportedOperationException("Not supported yet.");
	}

	/**
	 * Getter for property hessian.
	 * 
	 * @return Value of property hessian.
	 */
	@Override
	public Matrix getHessian() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}