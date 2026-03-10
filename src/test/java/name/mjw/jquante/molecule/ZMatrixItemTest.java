package name.mjw.jquante.molecule;

import static org.junit.jupiter.api.Assertions.*;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

class ZMatrixItemTest {

    @Test
    void constructorSetsReferenceAtomAndValue() {
        Atom atom = new Atom("H", new Vector3D(0, 0, 0));
        ZMatrixItem item = new ZMatrixItem(atom, 1.5);
        assertSame(atom, item.getReferenceAtom());
        assertEquals(1.5, item.getValue(), 1e-10);
    }

    @Test
    void setReferenceAtom() {
        Atom a1 = new Atom("H", new Vector3D(0, 0, 0));
        Atom a2 = new Atom("O", new Vector3D(1, 0, 0));
        ZMatrixItem item = new ZMatrixItem(a1, 1.0);
        item.setReferenceAtom(a2);
        assertSame(a2, item.getReferenceAtom());
    }

    @Test
    void setValue() {
        Atom atom = new Atom("H", new Vector3D(0, 0, 0));
        ZMatrixItem item = new ZMatrixItem(atom, 1.0);
        item.setValue(2.5);
        assertEquals(2.5, item.getValue(), 1e-10);
    }

    @Test
    void cloneCreatesIndependentCopyWithSameData() throws CloneNotSupportedException {
        Atom atom = new Atom("C", new Vector3D(1, 2, 3));
        ZMatrixItem original = new ZMatrixItem(atom, 3.14);
        ZMatrixItem clone = (ZMatrixItem) original.clone();

        assertNotSame(original, clone);
        assertEquals(original.getValue(), clone.getValue(), 1e-10);
        assertSame(original.getReferenceAtom(), clone.getReferenceAtom());
    }

    @Test
    void cloneValueIsIndependent() throws CloneNotSupportedException {
        Atom atom = new Atom("N", new Vector3D(0, 0, 0));
        ZMatrixItem original = new ZMatrixItem(atom, 1.0);
        ZMatrixItem clone = (ZMatrixItem) original.clone();

        clone.setValue(99.0);
        assertEquals(1.0, original.getValue(), 1e-10);
    }

    @Test
    void toStringContainsSymbolAndValue() {
        Atom atom = new Atom("H", new Vector3D(0, 0, 0));
        ZMatrixItem item = new ZMatrixItem(atom, 1.5);
        String s = item.toString();
        assertTrue(s.contains("H"));
        assertTrue(s.contains("1.5"));
    }
}
