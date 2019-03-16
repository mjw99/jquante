package name.mjw.jquante.math.qm.integral;

import static org.junit.jupiter.api.Assertions.*;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.BeforeAll;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.BasisSetLibrary;
import name.mjw.jquante.math.qm.OneElectronIntegrals;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.impl.MoleculeImpl;

class OneElectronIntegralsTest {

	double diff = 0.0001;
	static OneElectronIntegrals e1 = null;

	@BeforeAll
	public static void setUp() {

		// Create molecule
		Atom O = new Atom("O", 6.0, new Vector3D(0.00000000, 0.000000, 0.119748));
		Atom H1 = new Atom("H", 1.0, new Vector3D(0.00000000, 0.761561, -0.478993));
		Atom H2 = new Atom("H", 1.0, new Vector3D(0.00000000, -0.761561, -0.478993));

		Molecule water = new MoleculeImpl("water");
		water.addAtom(O);
		water.addAtom(H1);
		water.addAtom(H2);

		// Read Basis
		BasisSetLibrary bf = null;

		try {
			bf = new BasisSetLibrary(water, "sto-3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// compute integrals
		e1 = new OneElectronIntegrals(bf, water);

	}

	@Test
	void testHCore() {
		assertEquals(-32.70775590590571, e1.getHCore().getEntry(0, 0), diff);
	}

	@Test
	void testgetOverlap() {
		assertEquals(1.0, e1.getOverlap().getEntry(0, 0), diff);
	}

}
