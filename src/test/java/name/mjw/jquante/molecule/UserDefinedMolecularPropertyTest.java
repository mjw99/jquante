package name.mjw.jquante.molecule;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class UserDefinedMolecularPropertyTest {

    @Test
    void constructorSetsNameAndValue() {
        UserDefinedMolecularProperty p = new UserDefinedMolecularProperty("energy", -75.5);
        assertEquals("energy", p.getName());
        assertEquals(-75.5, p.getValue());
    }

    @Test
    void defaultConstructorHasNullFields() {
        UserDefinedMolecularProperty p = new UserDefinedMolecularProperty();
        assertNull(p.getName());
        assertNull(p.getValue());
    }

    @Test
    void setName() {
        UserDefinedMolecularProperty p = new UserDefinedMolecularProperty();
        p.setName("dipole");
        assertEquals("dipole", p.getName());
    }

    @Test
    void setValue() {
        UserDefinedMolecularProperty p = new UserDefinedMolecularProperty();
        p.setValue(3.14);
        assertEquals(3.14, p.getValue());
    }

    @Test
    void originDependentDefaultIsFalse() {
        UserDefinedMolecularProperty p = new UserDefinedMolecularProperty("x", 0);
        assertFalse(p.isOriginDependent());
    }

    @Test
    void setOriginDependent() {
        UserDefinedMolecularProperty p = new UserDefinedMolecularProperty("x", 0);
        p.setOriginDependent(true);
        assertTrue(p.isOriginDependent());
    }
}
