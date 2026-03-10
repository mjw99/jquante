package name.mjw.jquante.config.impl;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class AtomPropertyTest {

    @Test
    void constructorSetsName() {
        AtomProperty ap = new AtomProperty("Hydrogen", 1, 1.008);
        assertEquals("Hydrogen", ap.getName());
    }

    @Test
    void constructorSetsAtomicNumber() {
        AtomProperty ap = new AtomProperty("Hydrogen", 1, 1.008);
        assertEquals(1, ap.getAtomicNumber());
    }

    @Test
    void constructorSetsAtomicWeight() {
        AtomProperty ap = new AtomProperty("Hydrogen", 1, 1.008);
        assertEquals(1.008, ap.getAtomicWeight(), 1e-10);
    }

    @Test
    void setName() {
        AtomProperty ap = new AtomProperty("X", 0, 0.0);
        ap.setName("Carbon");
        assertEquals("Carbon", ap.getName());
    }

    @Test
    void setAtomicNumber() {
        AtomProperty ap = new AtomProperty("X", 0, 0.0);
        ap.setAtomicNumber(6);
        assertEquals(6, ap.getAtomicNumber());
    }

    @Test
    void setAtomicWeight() {
        AtomProperty ap = new AtomProperty("X", 0, 0.0);
        ap.setAtomicWeight(12.011);
        assertEquals(12.011, ap.getAtomicWeight(), 1e-10);
    }

    @Test
    void getValueReturnsSelf() {
        AtomProperty ap = new AtomProperty("H", 1, 1.008);
        assertSame(ap, ap.getValue());
    }

    @Test
    void setValueCopiesFromOther() {
        AtomProperty target = new AtomProperty("X", 0, 0.0);
        AtomProperty source = new AtomProperty("Oxygen", 8, 15.999);
        target.setValue(source);
        assertEquals("Oxygen", target.getName());
        assertEquals(8, target.getAtomicNumber());
        assertEquals(15.999, target.getAtomicWeight(), 1e-10);
    }
}
