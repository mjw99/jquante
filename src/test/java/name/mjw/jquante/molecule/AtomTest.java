package name.mjw.jquante.molecule;

import static org.junit.jupiter.api.Assertions.*;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

class AtomTest {

    private Atom atom;

    @BeforeEach
    void setUp() {
        atom = new Atom("C", new Vector3D(1.0, 2.0, 3.0));
    }

    @Test
    void symbolIsCapitalized() {
        assertEquals("C", atom.getSymbol());
    }

    @Test
    void symbolFromLowercaseIsCapitalized() {
        Atom a = new Atom("h", new Vector3D(0, 0, 0));
        assertEquals("H", a.getSymbol());
    }

    @Test
    void atomCenterCoordinates() {
        assertEquals(1.0, atom.getX(), 1e-10);
        assertEquals(2.0, atom.getY(), 1e-10);
        assertEquals(3.0, atom.getZ(), 1e-10);
    }

    @Test
    void defaultAtomCenterUnitsIsAngstrom() {
        assertEquals(name.mjw.jquante.common.Units.ANGSTROM, atom.getAtomCenterUnits());
    }

    @Test
    void distanceFromSamePoint() {
        assertEquals(0.0, atom.distanceFrom(new Vector3D(1.0, 2.0, 3.0)), 1e-10);
    }

    @Test
    void distanceFromOrigin() {
        Atom origin = new Atom("H", new Vector3D(0, 0, 0));
        Atom far = new Atom("H", new Vector3D(3, 4, 0));
        assertEquals(5.0, far.distanceFrom(origin), 1e-10);
    }

    @Test
    void distanceFromVector() {
        assertEquals(0.0, atom.distanceFrom(new Vector3D(1.0, 2.0, 3.0)), 1e-10);
    }

    @Test
    void defaultIndexIsZero() {
        assertEquals(0, atom.getIndex());
    }

    @Test
    void indexWithConstructor() {
        Atom a = new Atom("N", new Vector3D(0, 0, 0), 5);
        assertEquals(5, a.getIndex());
    }

    @Test
    void addAndQueryConnection() {
        atom.addConnection(1, BondType.SINGLE_BOND);
        assertTrue(atom.isConnected(1));
    }

    @Test
    void notConnectedByDefault() {
        assertFalse(atom.isConnected(99));
    }

    @Test
    void getConnectivityReturnsBondType() {
        atom.addConnection(2, BondType.DOUBLE_BOND);
        assertEquals(BondType.DOUBLE_BOND, atom.getConnectivity(2));
    }

    @Test
    void getConnectivityReturnsNoBondForMissingAtom() {
        assertEquals(BondType.NO_BOND, atom.getConnectivity(99));
    }

    @Test
    void removeConnection() {
        atom.addConnection(3, BondType.SINGLE_BOND);
        atom.removeConnection(3);
        assertFalse(atom.isConnected(3));
    }

    @Test
    void addNoBondRemovesExistingConnection() {
        atom.addConnection(4, BondType.SINGLE_BOND);
        atom.addConnection(4, BondType.NO_BOND);
        assertFalse(atom.isConnected(4));
    }

    @Test
    void removeAllConnections() {
        atom.addConnection(1, BondType.SINGLE_BOND);
        atom.addConnection(2, BondType.DOUBLE_BOND);
        atom.removeAllConnections();
        assertEquals(0, atom.getDegree());
    }

    @Test
    void getDegreeCountsConnections() {
        atom.addConnection(1, BondType.SINGLE_BOND);
        atom.addConnection(2, BondType.SINGLE_BOND);
        assertEquals(2, atom.getDegree());
    }

    @Test
    void numberOfStrongBonds() {
        atom.addConnection(1, BondType.SINGLE_BOND);
        atom.addConnection(2, BondType.DOUBLE_BOND);
        assertEquals(2, atom.getNumberOfStrongBonds());
    }

    @Test
    void numberOfDoubleBonds() {
        atom.addConnection(1, BondType.SINGLE_BOND);
        atom.addConnection(2, BondType.DOUBLE_BOND);
        assertEquals(1, atom.getNumberOfDoubleBonds());
    }

    @Test
    void hasResonantBondsWithAromaticBond() {
        atom.addConnection(1, BondType.AROMATIC_BOND);
        assertTrue(atom.hasResonantBonds());
    }

    @Test
    void noResonantBondsWithSingleBond() {
        atom.addConnection(1, BondType.SINGLE_BOND);
        assertFalse(atom.hasResonantBonds());
    }

    @Test
    void equalsSameAtom() {
        Atom a = new Atom("C", new Vector3D(1.0, 2.0, 3.0));
        assertEquals(a, atom);
    }

    @Test
    void notEqualsDifferentSymbol() {
        Atom a = new Atom("N", new Vector3D(1.0, 2.0, 3.0));
        assertNotEquals(a, atom);
    }

    @Test
    void hashCodeConsistentWithEquals() {
        Atom a = new Atom("C", new Vector3D(1.0, 2.0, 3.0));
        assertEquals(a.hashCode(), atom.hashCode());
    }

    @Test
    void toStringContainsSymbol() {
        assertTrue(atom.toString().contains("C"));
    }
}
