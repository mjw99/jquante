package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.Iterables;
import com.google.common.collect.Multimap;

import name.mjw.jquante.math.qm.basis.AtomicBasis;
import name.mjw.jquante.math.qm.basis.BasisReader;
import name.mjw.jquante.math.qm.basis.BasisSet;
import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Orbital;
import name.mjw.jquante.math.qm.basis.Power;
import name.mjw.jquante.math.qm.basis.PowerList;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.UserDefinedAtomProperty;
import name.mjw.jquante.molecule.event.MoleculeStateChangeEvent;
import name.mjw.jquante.molecule.event.MoleculeStateChangeListener;

/**
 * Class to construct basis functions of a given molecule and a basis set
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class BasisFunctions {

	private static final Logger LOG = LogManager.getLogger(BasisFunctions.class);

	/**
	 * Holds value of property basisFunctions.
	 */
	private ArrayList<ContractedGaussian> basisFunctions;

	private Multimap<Integer, ContractedGaussian> shells;

	private Multimap<ContractedGaussian, ContractedGaussian> shellPairs;

	private Molecule molecule;

	private MoleculeStateChangeListener molStateChangeListener;

	/**
	 * Creates a new instance of BasisFunctions.
	 * 
	 * @param molecule  the Molecule whose basis function is requested.
	 * @param basisName the name of the basis set (like sto3g).
	 * @throws Exception the basisName was not found.
	 *
	 */
	public BasisFunctions(Molecule molecule, String basisName) throws Exception {
		// initialise the basis functions
		getBasisFunctions(molecule, basisName);
		this.basisName = basisName;
		this.molecule = molecule;

		// and initialise the shell list
		initShellList();

		// Shell Pair list
		initShellPairList();

		molStateChangeListener = new MoleculeStateChangeListener() {
			@Override
			public void moleculeChanged(MoleculeStateChangeEvent event) {
				try {
					getBasisFunctions(BasisFunctions.this.molecule, BasisFunctions.this.basisName);
					initShellList();
				} catch (Exception e) {
					LOG.error("Unable to update basis function! ");
				}
			}
		};
		molecule.addMoleculeStateChangeListener(molStateChangeListener);
	}

	private String basisName;

	/**
	 * Get the value of basisName
	 * 
	 * @return the value of basisName
	 */
	public String getBasisName() {
		return basisName;
	}

	/**
	 * Set the value of basisName
	 * 
	 * @param basisName new value of basisName
	 */
	public void setBasisName(String basisName) {
		this.basisName = basisName;
	}

	/**
	 * Getter for property basisFunctions.
	 * 
	 * @return Value of property basisFunctions.
	 */
	public List<ContractedGaussian> getBasisFunctions() {
		return this.basisFunctions;
	}

	public Multimap<Integer, ContractedGaussian> getShells() {
		return shells;
	}

	public Multimap<ContractedGaussian, ContractedGaussian> getShellPairs() {
		return shellPairs;
	}

	/**
	 * Getter for property basisFunctions.
	 * 
	 * @param molecule  the Molecule whose basis function is requested
	 * @param basisName the name of the basis set (like sto3g)
	 * @return Value of property basisFunctions.
	 */
	private ArrayList<ContractedGaussian> getBasisFunctions(Molecule molecule, String basisName) throws Exception {
		BasisSet basis = BasisReader.getInstance().readBasis(basisName);
		Iterator<Atom> atoms = molecule.getAtoms();

		basisFunctions = new ArrayList<>();

		Atom atom;
		AtomicBasis atomicBasis;
		while (atoms.hasNext()) { // loop over atoms
			atom = atoms.next();
			atomicBasis = basis.getAtomicBasis(atom.getSymbol());

			Iterator<Orbital> orbitals = atomicBasis.getOrbitals().iterator();
			Orbital orbital;
			ArrayList<ContractedGaussian> atomicFunctions = new ArrayList<>();

			while (orbitals.hasNext()) { // loop over atom orbitals
				orbital = orbitals.next();

				Iterator<Power> pList = PowerList.getInstance().getPowerList(orbital.getType());
				Power power;
				while (pList.hasNext()) { // and the power list, sp2 etc..
					power = pList.next();

					ContractedGaussian cg = new ContractedGaussian(atom, power);
					Iterator<Double> coeff = orbital.getCoefficients().iterator();
					Iterator<Double> exp = orbital.getExponents().iterator();

					while (coeff.hasNext()) { // build the CG from PGs
						cg.addPrimitive(exp.next().doubleValue(), coeff.next().doubleValue());
					} // end while

					cg.normalize();
					cg.setIndex(basisFunctions.size()); // send an index
					basisFunctions.add(cg); // add this CG to list
					atomicFunctions.add(cg); // add the reference to atom list
				}
			}

			// save a reference of the basis functions centered on
			// this atom as a user defined property of the atom
			try {
				atom.addUserDefinedAtomProperty(new UserDefinedAtomProperty("basisFunctions", atomicFunctions));
			} catch (UnsupportedOperationException e) {
				UserDefinedAtomProperty up = atom.getUserDefinedAtomProperty("basisFunctions");
				up.setValue(atomicFunctions);
			}
		}

		return this.basisFunctions;
	}

	/**
	 * Initialise the shell list
	 */
	private void initShellList() {
		shells = ArrayListMultimap.create();

		for (int i = 0; i < basisFunctions.size(); i++) {

			Boolean addToExistingShell = false;
			Integer existingShellIndex = 0;

			// First entry always goes in.
			if (i == 0) {

				shells.put(0, basisFunctions.get(0));

			} else {

				// Look to see if the basisfunction is already present in the known shells
				Iterator<Integer> keyIterator = shells.keys().iterator();

				while (keyIterator.hasNext()) {
					Integer key = keyIterator.next();
					Collection<ContractedGaussian> shellCgs = shells.get(key);

					for (ContractedGaussian shellCg : shellCgs) {
						if (basisFunctions.get(i).isSameShell(shellCg)) {
							// Flag this to the existing shell.
							// This is a workaround for Multimap checkForComodification exceptions.
							addToExistingShell = true;
							existingShellIndex = key;
						}
					}

				}

				if (addToExistingShell) {
					shells.put(existingShellIndex, basisFunctions.get(i));
				} else {
					int max = shells.keys().stream().mapToInt(v -> v).max().getAsInt();
					shells.put(max + 1, basisFunctions.get(i));
				}

			}

		}

	}

	private void initShellPairList() {

		shellPairs = ArrayListMultimap.create();

		for (int i = 0; i < shells.keySet().size(); i++) {
			for (int j = 0; j <= i; j++) {
				shellPairs.put(Iterables.get(shells.get(i), 0), Iterables.get(shells.get(j), 0));
			}
		}

	}

}
