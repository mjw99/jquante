package name.mjw.jquante.math.qm;

import static org.junit.Assert.*;

import name.mjw.jquante.math.geom.Point3D;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.impl.MoleculeImpl;

import org.apache.log4j.Logger;
import org.junit.Test;

public class SCFMethodFactoryTest {

	private final Logger LOG = Logger.getLogger(SCFMethodFactoryTest.class);

	double diff = 0.00001;

	@Test
	public void SinglePointHFWaterSTO3G() {

		// Create molecule
		Atom O = new Atom("O", 6.0, new Point3D(0.00000000, 0.000000, 0.119748));
		Atom H1 = new Atom("H", 1.0, new Point3D(0.00000000, 0.761561,
				-0.478993));
		Atom H2 = new Atom("H", 1.0, new Point3D(0.00000000, -0.761561,
				-0.478993));

		Molecule water = new MoleculeImpl("water");
		water.addAtom(H1);
		water.addAtom(O);
		water.addAtom(H2);

		long t1 = System.currentTimeMillis();
		// Read Basis
		BasisFunctions bf = null;

		try {
			bf = new BasisFunctions(water, "sto3g");

		} catch (Exception e) {

			e.printStackTrace();
		}

		// compute integrals
		OneElectronIntegrals e1 = new OneElectronIntegrals(bf, water);
		TwoElectronIntegrals e2 = new TwoElectronIntegrals(bf);

		long t2 = System.currentTimeMillis();

		// do SCF
		SCFMethod scfm = SCFMethodFactory.getInstance().getSCFMethod(water, e1,
				e2, SCFType.HARTREE_FOCK);
		scfm.scf();

		long t3 = System.currentTimeMillis();

		LOG.debug("Time till 2E : " + (t2 - t1) + " ms");
		LOG.debug("Time for SCF : " + (t3 - t2) + " ms");

		assertEquals(-74.9645174916, scfm.getEnergy(), diff);

		// orbital energies
		double[] ev = scfm.getOrbE();

		assertEquals(-20.244508003283794, ev[0], diff);
		assertEquals(-1.2637569211246862, ev[1], diff);
		assertEquals(-0.6106330109275634, ev[2], diff);
		assertEquals(-0.4535337774289125, ev[3], diff);
		assertEquals(-0.3913211137352477, ev[4], diff);
		assertEquals(0.5958982943824441, ev[5], diff);
		assertEquals(0.7260124013928394, ev[6], diff);

	}
}
