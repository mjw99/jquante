package name.mjw.jquante.molecule;

import static org.junit.jupiter.api.Assertions.*;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.molecule.impl.MoleculeImpl;
import name.mjw.jquante.test.Fixtures;

class MolecularFormulaTest {

    @Test
    void waterFormulaHasTwoHydrogen() {
        MolecularFormula formula = new MolecularFormula(Fixtures.getWater());
        assertEquals(2, formula.getAtomCount("H"));
    }

    @Test
    void waterFormulaHasOneOxygen() {
        MolecularFormula formula = new MolecularFormula(Fixtures.getWater());
        assertEquals(1, formula.getAtomCount("O"));
    }

    @Test
    void waterFormulaString() {
        MolecularFormula formula = new MolecularFormula(Fixtures.getWater());
        // elements are sorted alphabetically: H before O
        assertTrue(formula.getFormula().contains("H2"));
        assertTrue(formula.getFormula().contains("O1"));
    }

    @Test
    void hydrogenFormulaHasTwoHydrogen() {
        MolecularFormula formula = new MolecularFormula(Fixtures.getHydrogen());
        assertEquals(2, formula.getAtomCount("H"));
    }

    @Test
    void atomCountZeroForMissingElement() {
        MolecularFormula formula = new MolecularFormula(Fixtures.getWater());
        assertEquals(0, formula.getAtomCount("N"));
    }

    @Test
    void toStringMatchesGetFormula() {
        MolecularFormula formula = new MolecularFormula(Fixtures.getWater());
        assertEquals(formula.getFormula(), formula.toString());
    }

    @Test
    void addFormulasProducesCorrectCounts() {
        MolecularFormula f1 = new MolecularFormula(Fixtures.getWater());
        MolecularFormula f2 = new MolecularFormula(Fixtures.getWater());
        MolecularFormula sum = f1.add(f2);
        assertEquals(4, sum.getAtomCount("H"));
        assertEquals(2, sum.getAtomCount("O"));
    }

    @Test
    void subtractFormulasProducesCorrectCounts() {
        MolecularFormula f1 = new MolecularFormula(Fixtures.getWater());
        MolecularFormula f2 = new MolecularFormula(Fixtures.getHydrogen());
        MolecularFormula diff = f1.sub(f2);
        assertEquals(0, diff.getAtomCount("H"));
        assertEquals(1, diff.getAtomCount("O"));
    }

    @Test
    void equalFormulasAreEqual() {
        MolecularFormula f1 = new MolecularFormula(Fixtures.getWater());
        MolecularFormula f2 = new MolecularFormula(Fixtures.getWater());
        assertEquals(f1, f2);
    }

    @Test
    void differentFormulasAreNotEqual() {
        MolecularFormula f1 = new MolecularFormula(Fixtures.getWater());
        MolecularFormula f2 = new MolecularFormula(Fixtures.getHydrogen());
        assertNotEquals(f1, f2);
    }

    @Test
    void formulaUpdatesWhenAtomAdded() {
        MoleculeImpl mol = new MoleculeImpl("test");
        mol.addAtom("H", new Vector3D(0, 0, 0));
        MolecularFormula formula = new MolecularFormula(mol);
        assertEquals(1, formula.getAtomCount("H"));

        mol.addAtom("H", new Vector3D(1, 0, 0));
        assertEquals(2, formula.getAtomCount("H"));
    }
}
