package name.mjw.jquante.test;

import org.hipparchus.geometry.euclidean.threed.Vector3D;

import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.impl.MoleculeImpl;

public class Fixtures {

	static public Molecule getWaterMolecule() {
		Atom O = new Atom("O", new Vector3D(0.0000000, 0.000000, 0.119748));
		Atom H1 = new Atom("H", new Vector3D(0.00000000, 0.761561, -0.478993));
		Atom H2 = new Atom("H", new Vector3D(0.00000000, -0.761561, -0.478993));

		Molecule water = new MoleculeImpl("water");
		water.addAtom(O);
		water.addAtom(H1);
		water.addAtom(H2);

		return water;
	}

	static public Molecule getHydrogenMolecule() {

		Atom H1 = new Atom("H", new Vector3D(0.00000000, 0.00000000, 0.00000000));
		Atom H2 = new Atom("H", new Vector3D(0.74000000, 0.00000000, 0.00000000));

		Molecule hydrogen = new MoleculeImpl("hydrogen");
		hydrogen.addAtom(H1);
		hydrogen.addAtom(H2);

		return hydrogen;
	}

	static public Molecule getHydrogenFluoride() {
		// Create molecule
		Atom H = new Atom("H", new Vector3D(0.00000000, 0.00000000, 0.00000000));
		Atom F = new Atom("F", new Vector3D(0.91700000, 0.00000000, 0.00000000));

		Molecule hydrogenFluoride = new MoleculeImpl("hydrogenFluoride");
		hydrogenFluoride.addAtom(H);
		hydrogenFluoride.addAtom(F);

		return hydrogenFluoride;

	}

}
