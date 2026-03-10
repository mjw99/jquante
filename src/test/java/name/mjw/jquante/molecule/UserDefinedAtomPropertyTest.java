package name.mjw.jquante.molecule;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class UserDefinedAtomPropertyTest {

    @Test
    void constructorSetsNameAndValue() {
        UserDefinedAtomProperty p = new UserDefinedAtomProperty("charge", 1.5);
        assertEquals("charge", p.getName());
        assertEquals(1.5, p.getValue());
    }

    @Test
    void defaultConstructorHasNullFields() {
        UserDefinedAtomProperty p = new UserDefinedAtomProperty();
        assertNull(p.getName());
        assertNull(p.getValue());
    }

    @Test
    void setName() {
        UserDefinedAtomProperty p = new UserDefinedAtomProperty();
        p.setName("energy");
        assertEquals("energy", p.getName());
    }

    @Test
    void setValue() {
        UserDefinedAtomProperty p = new UserDefinedAtomProperty();
        p.setValue("some value");
        assertEquals("some value", p.getValue());
    }

    @Test
    void valueCanBeAnySerializable() {
        UserDefinedAtomProperty p = new UserDefinedAtomProperty("count", 42);
        assertEquals(42, p.getValue());
    }
}
