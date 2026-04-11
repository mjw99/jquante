package name.mjw.jquante.molecule;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class BondTypeTest {

    private static final double DIFF = 1e-10;

    @Test
    void singleBondOrder() {
        assertEquals(1.0, BondType.SINGLE_BOND.getBondOrder(), DIFF);
    }

    @Test
    void doubleBondOrder() {
        assertEquals(2.0, BondType.DOUBLE_BOND.getBondOrder(), DIFF);
    }

    @Test
    void tripleBondOrder() {
        assertEquals(3.0, BondType.TRIPLE_BOND.getBondOrder(), DIFF);
    }

    @Test
    void quadrupleBondOrder() {
        assertEquals(4.0, BondType.QUADRUPLE_BOND.getBondOrder(), DIFF);
    }

    @Test
    void aromaticBondOrder() {
        assertEquals(1.5, BondType.AROMATIC_BOND.getBondOrder(), DIFF);
    }

    @Test
    void amideBondOrder() {
        assertEquals(1.41, BondType.AMIDE_BOND.getBondOrder(), DIFF);
    }

    @Test
    void noBondOrderIsZero() {
        assertEquals(0.0, BondType.NO_BOND.getBondOrder(), DIFF);
    }

    @Test
    void weakBondOrderIsZero() {
        assertEquals(0.0, BondType.WEAK_BOND.getBondOrder(), DIFF);
    }

    @Test
    void ionicBondOrderIsZero() {
        assertEquals(0.0, BondType.IONIC_BOND.getBondOrder(), DIFF);
    }

    @Test
    void getBondTypeForSingleBond() {
        assertEquals(BondType.SINGLE_BOND, BondType.getBondTypeFor("Single bond"));
    }

    @Test
    void getBondTypeForDoubleBond() {
        assertEquals(BondType.DOUBLE_BOND, BondType.getBondTypeFor("Double bond"));
    }

    @Test
    void getBondTypeForTripleBond() {
        assertEquals(BondType.TRIPLE_BOND, BondType.getBondTypeFor("Triple bond"));
    }

    @Test
    void getBondTypeForAromaticBond() {
        assertEquals(BondType.AROMATIC_BOND, BondType.getBondTypeFor("Aromatic bond"));
    }

    @Test
    void getBondTypeForUnknownReturnsNoBond() {
        assertEquals(BondType.NO_BOND, BondType.getBondTypeFor("Unknown type"));
    }

    @Test
    void getBondTypeForEmptyStringReturnsNoBond() {
        assertEquals(BondType.NO_BOND, BondType.getBondTypeFor(""));
    }

    @Test
    void isStrongBondTrueForSingleBond() {
        assertTrue(BondType.SINGLE_BOND.isStrongBond());
    }

    @Test
    void isStrongBondTrueForDoubleBond() {
        assertTrue(BondType.DOUBLE_BOND.isStrongBond());
    }

    @Test
    void isStrongBondTrueForTripleBond() {
        assertTrue(BondType.TRIPLE_BOND.isStrongBond());
    }

    @Test
    void isStrongBondTrueForAromaticBond() {
        assertTrue(BondType.AROMATIC_BOND.isStrongBond());
    }

    @Test
    void isStrongBondTrueForAmideBond() {
        assertTrue(BondType.AMIDE_BOND.isStrongBond());
    }

    @Test
    void isStrongBondFalseForNoBond() {
        assertFalse(BondType.NO_BOND.isStrongBond());
    }

    @Test
    void isStrongBondFalseForWeakBond() {
        assertFalse(BondType.WEAK_BOND.isStrongBond());
    }

    @Test
    void toStringForSingleBond() {
        assertEquals("Single bond", BondType.SINGLE_BOND.toString());
    }

    @Test
    void toStringForDoubleBond() {
        assertEquals("Double bond", BondType.DOUBLE_BOND.toString());
    }

    @Test
    void toStringForNoBond() {
        assertEquals("No bond", BondType.NO_BOND.toString());
    }

    @Test
    void toStringForIonicBond() {
        assertEquals("Ionic bond", BondType.IONIC_BOND.toString());
    }
}
