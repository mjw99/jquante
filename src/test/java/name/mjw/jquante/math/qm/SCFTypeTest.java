package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class SCFTypeTest {

    @Test
    void hartreeFockTypeValue() {
        assertEquals(1, SCFType.HARTREE_FOCK.getType());
    }

    @Test
    void hartreeFockDirectTypeValue() {
        assertEquals(2, SCFType.HARTREE_FOCK_DIRECT.getType());
    }

    @Test
    void mollerPlessetTypeValue() {
        assertEquals(3, SCFType.MOLLER_PLESSET.getType());
    }

    @Test
    void equalsReflexive() {
        assertTrue(SCFType.HARTREE_FOCK.equals(SCFType.HARTREE_FOCK));
    }

    @Test
    void equalsSameInstance() {
        SCFType hf = SCFType.HARTREE_FOCK;
        assertTrue(hf.equals(hf));
    }

    @Test
    void notEqualsForDifferentTypes() {
        assertFalse(SCFType.HARTREE_FOCK.equals(SCFType.HARTREE_FOCK_DIRECT));
        assertFalse(SCFType.HARTREE_FOCK_DIRECT.equals(SCFType.HARTREE_FOCK));
    }

    @Test
    void notEqualsForString() {
        assertFalse(SCFType.HARTREE_FOCK.equals("HARTREE_FOCK"));
    }

    @Test
    void notEqualsForNull() {
        assertFalse(SCFType.HARTREE_FOCK.equals(null));
    }

    @Test
    void hashCodeIsConsistent() {
        int h1 = SCFType.HARTREE_FOCK.hashCode();
        int h2 = SCFType.HARTREE_FOCK.hashCode();
        assertEquals(h1, h2);
    }

    @Test
    void equalObjectsHaveSameHashCode() {
        // Two SCFType objects with the same type value should be equal and have same hashCode
        assertEquals(SCFType.MOLLER_PLESSET.hashCode(), SCFType.UNRESTRICTED_HARTREE_FOCK_DIRECT.hashCode());
    }

    @Test
    void toStringForHartreeFockContainsHartreeFock() {
        String s = SCFType.HARTREE_FOCK.toString();
        assertTrue(s.contains("Hartree Fock"), "Expected 'Hartree Fock' in: " + s);
    }

    @Test
    void toStringForHartreeFockDirectContainsDirect() {
        String s = SCFType.HARTREE_FOCK_DIRECT.toString();
        assertTrue(s.contains("direct") || s.contains("Direct"),
                "Expected 'direct' or 'Direct' in: " + s);
    }

    @Test
    void toStringForHartreeFockDoesNotContainDirect() {
        String s = SCFType.HARTREE_FOCK.toString();
        assertFalse(s.contains("direct") || s.contains("Direct"),
                "HARTREE_FOCK toString should not contain 'direct'");
    }
}
