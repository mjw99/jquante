package name.mjw.jquante.molecule.impl;

import static org.junit.jupiter.api.Assertions.*;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.geom.BoundingBox;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.BondType;
import name.mjw.jquante.test.Fixtures;

class MoleculeImplTest {

    private MoleculeImpl molecule;

    @BeforeEach
    void setUp() {
        molecule = new MoleculeImpl("test");
    }

    @Test
    void titleIsSet() {
        assertEquals("test", molecule.getTitle());
    }

    @Test
    void defaultTitleIsMolecule() {
        assertEquals("Molecule", new MoleculeImpl().getTitle());
    }

    @Test
    void emptyMoleculeHasZeroAtoms() {
        assertEquals(0, molecule.getNumberOfAtoms());
    }

    @Test
    void addAtomIncreasesCount() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        assertEquals(1, molecule.getNumberOfAtoms());
    }

    @Test
    void addAtomAssignsIndex() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        molecule.addAtom("O", new Vector3D(1, 0, 0));
        assertEquals(0, molecule.getAtom(0).getIndex());
        assertEquals(1, molecule.getAtom(1).getIndex());
    }

    @Test
    void addAtomByCoordinates() {
        molecule.addAtom("C", 1.0, 2.0, 3.0);
        Atom atom = molecule.getAtom(0);
        assertEquals(1.0, atom.getX(), 1e-10);
        assertEquals(2.0, atom.getY(), 1e-10);
        assertEquals(3.0, atom.getZ(), 1e-10);
    }

    @Test
    void removeAtomDecreasesCount() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        molecule.addAtom("H", new Vector3D(1, 0, 0));
        molecule.removeAtomAt(0);
        assertEquals(1, molecule.getNumberOfAtoms());
    }

    @Test
    void removeAtomReindexes() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        molecule.addAtom("O", new Vector3D(1, 0, 0));
        molecule.removeAtomAt(0);
        assertEquals(0, molecule.getAtom(0).getIndex());
    }

    @Test
    void numberOfElectronsForWater() {
        // O=8, H=1, H=1 => 10 electrons
        assertEquals(10, Fixtures.getWater().getNumberOfElectrons());
    }

    @Test
    void molecularMassForHydrogen() {
        // H2: 2 * ~1.008 amu
        double mass = Fixtures.getHydrogen().getMolecularMass();
        assertEquals(2.016, mass, 0.01);
    }

    @Test
    void setBondTypeAndQuery() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        molecule.addAtom("H", new Vector3D(1, 0, 0));
        molecule.setBondType(0, 1, BondType.SINGLE_BOND);
        assertTrue(molecule.isBonded(0, 1));
        assertTrue(molecule.isBonded(1, 0));
    }

    @Test
    void removeBondBetween() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        molecule.addAtom("H", new Vector3D(1, 0, 0));
        molecule.setBondType(0, 1, BondType.SINGLE_BOND);
        molecule.removeBondBetween(0, 1);
        assertFalse(molecule.isBonded(0, 1));
    }

    @Test
    void getBondType() {
        molecule.addAtom("C", new Vector3D(0, 0, 0));
        molecule.addAtom("C", new Vector3D(1, 0, 0));
        molecule.setBondType(0, 1, BondType.DOUBLE_BOND);
        assertEquals(BondType.DOUBLE_BOND, molecule.getBondType(0, 1));
    }

    @Test
    void boundingBoxForWater() {
        BoundingBox bb = Fixtures.getWater().getBoundingBox();
        assertNotNull(bb);
        assertTrue(bb.getXWidth() >= 0.0);
        assertTrue(bb.getYWidth() > 0.0);
    }

    @Test
    void emptyMoleculeReturnsZeroBoundingBox() {
        BoundingBox bb = molecule.getBoundingBox();
        assertEquals(0.0, bb.getXWidth(), 1e-10);
    }

    @Test
    void formulaForWater() {
        assertEquals(2, Fixtures.getWater().getFormula().getAtomCount("H"));
        assertEquals(1, Fixtures.getWater().getFormula().getAtomCount("O"));
    }

    @Test
    void setTitleUpdatesTitle() {
        molecule.setTitle("ethanol");
        assertEquals("ethanol", molecule.getTitle());
    }

    @Test
    void soundsLikeFindsAtom() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        Atom query = new Atom("H", new Vector3D(0, 0, 0));
        assertEquals(0, molecule.soundsLike(query));
    }

    @Test
    void soundsLikeReturnsMinusOneForMissing() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        Atom query = new Atom("H", new Vector3D(10, 10, 10));
        assertEquals(-1, molecule.soundsLike(query));
    }

    @Test
    void resetAtomCoordinatesUpdatesPositions() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        molecule.addAtom("H", new Vector3D(1, 0, 0));
        double[] newCoords = {2.0, 0.0, 0.0, 3.0, 0.0, 0.0};
        molecule.resetAtomCoordinates(newCoords, false);
        assertEquals(2.0, molecule.getAtom(0).getX(), 1e-10);
        assertEquals(3.0, molecule.getAtom(1).getX(), 1e-10);
    }

    @Test
    void resetAtomCoordinatesThrowsForWrongSize() {
        molecule.addAtom("H", new Vector3D(0, 0, 0));
        assertThrows(IllegalArgumentException.class,
                () -> molecule.resetAtomCoordinates(new double[]{1.0, 2.0}, false));
    }
}
