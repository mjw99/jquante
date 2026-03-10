package name.mjw.jquante.config.impl;

import static org.junit.jupiter.api.Assertions.*;

import org.junit.jupiter.api.Test;

class AtomInfoTest {

    @Test
    void singletonReturnsSameInstance() {
        assertSame(AtomInfo.getInstance(), AtomInfo.getInstance());
    }

    @Test
    void hydrogenAtomicNumber() {
        assertEquals(1, AtomInfo.getInstance().getAtomicNumber("H"));
    }

    @Test
    void carbonAtomicNumber() {
        assertEquals(6, AtomInfo.getInstance().getAtomicNumber("C"));
    }

    @Test
    void nitrogenAtomicNumber() {
        assertEquals(7, AtomInfo.getInstance().getAtomicNumber("N"));
    }

    @Test
    void oxygenAtomicNumber() {
        assertEquals(8, AtomInfo.getInstance().getAtomicNumber("O"));
    }

    @Test
    void hydrogenAtomicWeight() {
        assertEquals(1.008, AtomInfo.getInstance().getAtomicWeight("H"), 0.01);
    }

    @Test
    void carbonAtomicWeight() {
        assertEquals(12.011, AtomInfo.getInstance().getAtomicWeight("C"), 0.01);
    }

    @Test
    void hydrogenName() {
        assertEquals("Hydrogen", AtomInfo.getInstance().getName("H"));
    }

    @Test
    void carbonName() {
        assertEquals("Carbon", AtomInfo.getInstance().getName("C"));
    }

    @Test
    void getSymbolByAtomicNumber() {
        assertEquals("H", AtomInfo.getInstance().getSymbol(1));
    }

    @Test
    void getSymbolByAtomicNumberCarbon() {
        assertEquals("C", AtomInfo.getInstance().getSymbol(6));
    }

    @Test
    void unknownSymbolReturnsX() {
        assertEquals(0, AtomInfo.getInstance().getAtomicNumber("Unobtainium"));
    }

    @Test
    void getParameterReturnsAtomProperty() {
        AtomProperty ap = (AtomProperty) AtomInfo.getInstance().getParameter("H");
        assertNotNull(ap);
        assertEquals(1, ap.getAtomicNumber());
    }
}
