package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.impl.MoleculeImpl;

/**
 * Validates the analytic gradient in HartreeFockForce against a central
 * finite-difference (numerical) gradient.  This test is independent of any
 * external reference, so it isolates whether the analytic implementation is
 * self-consistent.
 *
 * Displacement h = 0.001 Bohr (≈ 0.000529 Å).  We work in Angstroms for the
 * geometry and convert: h_angstrom = h_bohr * 0.52917724924.
 */
class HartreeFockForceNumericalTest {

    private static final double BOHR_TO_ANG = 0.52917724924;
    private static final double H_BOHR = 0.001;        // displacement in Bohr
    private static final double H_ANG  = H_BOHR * BOHR_TO_ANG;

    private double scfEnergy(Molecule mol) {
        BasisSetLibrary bsl;
        try {
            bsl = new BasisSetLibrary(mol, "sto-3g");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        OneElectronIntegrals e1 = new OneElectronIntegrals(bsl, mol);
        TwoElectronIntegrals e2 = new TwoElectronIntegrals(bsl);
        SCFMethod scfm = SCFMethodFactory.getInstance()
                .getSCFMethod(mol, e1, e2, SCFType.HARTREE_FOCK);
        scfm.scf();
        return scfm.getEnergy();
    }

    private SCFMethod runSCF(Molecule mol) {
        BasisSetLibrary bsl;
        try {
            bsl = new BasisSetLibrary(mol, "sto-3g");
        } catch (Exception e) {
            throw new RuntimeException(e);
        }
        OneElectronIntegrals e1 = new OneElectronIntegrals(bsl, mol);
        TwoElectronIntegrals e2 = new TwoElectronIntegrals(bsl);
        SCFMethod scfm = SCFMethodFactory.getInstance()
                .getSCFMethod(mol, e1, e2, SCFType.HARTREE_FOCK);
        scfm.scf();
        return scfm;
    }

    // -----------------------------------------------------------------------
    // H2: analytic force vs numerical finite-difference
    // H1 at (0, 0, 0), H2 at (0.74, 0, 0) Å
    // -----------------------------------------------------------------------

    @Test
    void h2H1AnalyticMatchesNumerical() {
        // Central difference for atom 0 (H1), x-displacement
        Molecule molP = new MoleculeImpl("h2p");
        molP.addAtom(new Atom("H", new Vector3D( H_ANG, 0, 0)));
        molP.addAtom(new Atom("H", new Vector3D(0.74,   0, 0)));

        Molecule molM = new MoleculeImpl("h2m");
        molM.addAtom(new Atom("H", new Vector3D(-H_ANG, 0, 0)));
        molM.addAtom(new Atom("H", new Vector3D(0.74,   0, 0)));

        // numerical gradient dE/dx (force = -dE/dx)
        double numGrad = (scfEnergy(molP) - scfEnergy(molM)) / (2.0 * H_BOHR);
        double numForce = -numGrad;

        // analytic force from jquante
        Molecule mol0 = new MoleculeImpl("h2");
        mol0.addAtom(new Atom("H", new Vector3D(0.0, 0, 0)));
        mol0.addAtom(new Atom("H", new Vector3D(0.74, 0, 0)));
        SCFMethod scfm = runSCF(mol0);
        double analyticForce = new HartreeFockForce().computeForce(0, scfm).getX();

        System.out.printf("H2 H1 force: numeric=%.6f  analytic=%.6f%n", numForce, analyticForce);
        assertEquals(numForce, analyticForce, 1e-3);
    }

    @Test
    void h2H2AnalyticMatchesNumerical() {
        Molecule molP = new MoleculeImpl("h2p");
        molP.addAtom(new Atom("H", new Vector3D(0,          0, 0)));
        molP.addAtom(new Atom("H", new Vector3D(0.74+H_ANG, 0, 0)));

        Molecule molM = new MoleculeImpl("h2m");
        molM.addAtom(new Atom("H", new Vector3D(0,          0, 0)));
        molM.addAtom(new Atom("H", new Vector3D(0.74-H_ANG, 0, 0)));

        double numGrad  = (scfEnergy(molP) - scfEnergy(molM)) / (2.0 * H_BOHR);
        double numForce = -numGrad;

        Molecule mol0 = new MoleculeImpl("h2");
        mol0.addAtom(new Atom("H", new Vector3D(0.0,  0, 0)));
        mol0.addAtom(new Atom("H", new Vector3D(0.74, 0, 0)));
        SCFMethod scfm = runSCF(mol0);
        double analyticForce = new HartreeFockForce().computeForce(1, scfm).getX();

        System.out.printf("H2 H2 force: numeric=%.6f  analytic=%.6f%n", numForce, analyticForce);
        assertEquals(numForce, analyticForce, 1e-3);
    }

    // -----------------------------------------------------------------------
    // Water: analytic force on O (z), H1 (y, z)
    // O at (0,0,0.119748), H1 at (0, 0.761561,-0.478993) Å
    // -----------------------------------------------------------------------

    @Test
    void waterOzAnalyticMatchesNumerical() {
        Molecule molP = new MoleculeImpl("waterP");
        molP.addAtom(new Atom("O", new Vector3D(0, 0,  0.119748 + H_ANG)));
        molP.addAtom(new Atom("H", new Vector3D(0,  0.761561, -0.478993)));
        molP.addAtom(new Atom("H", new Vector3D(0, -0.761561, -0.478993)));

        Molecule molM = new MoleculeImpl("waterM");
        molM.addAtom(new Atom("O", new Vector3D(0, 0,  0.119748 - H_ANG)));
        molM.addAtom(new Atom("H", new Vector3D(0,  0.761561, -0.478993)));
        molM.addAtom(new Atom("H", new Vector3D(0, -0.761561, -0.478993)));

        double numGrad  = (scfEnergy(molP) - scfEnergy(molM)) / (2.0 * H_BOHR);
        double numForce = -numGrad;

        Molecule mol0 = new MoleculeImpl("water");
        mol0.addAtom(new Atom("O", new Vector3D(0, 0, 0.119748)));
        mol0.addAtom(new Atom("H", new Vector3D(0,  0.761561, -0.478993)));
        mol0.addAtom(new Atom("H", new Vector3D(0, -0.761561, -0.478993)));
        SCFMethod scfm = runSCF(mol0);
        double analyticForce = new HartreeFockForce().computeForce(0, scfm).getZ();

        System.out.printf("Water O z force: numeric=%.6f  analytic=%.6f%n", numForce, analyticForce);
        assertEquals(numForce, analyticForce, 1e-3);
    }

    @Test
    void waterH1yAnalyticMatchesNumerical() {
        Molecule molP = new MoleculeImpl("waterP");
        molP.addAtom(new Atom("O", new Vector3D(0, 0, 0.119748)));
        molP.addAtom(new Atom("H", new Vector3D(0, 0.761561 + H_ANG, -0.478993)));
        molP.addAtom(new Atom("H", new Vector3D(0, -0.761561, -0.478993)));

        Molecule molM = new MoleculeImpl("waterM");
        molM.addAtom(new Atom("O", new Vector3D(0, 0, 0.119748)));
        molM.addAtom(new Atom("H", new Vector3D(0, 0.761561 - H_ANG, -0.478993)));
        molM.addAtom(new Atom("H", new Vector3D(0, -0.761561, -0.478993)));

        double numGrad  = (scfEnergy(molP) - scfEnergy(molM)) / (2.0 * H_BOHR);
        double numForce = -numGrad;

        Molecule mol0 = new MoleculeImpl("water");
        mol0.addAtom(new Atom("O", new Vector3D(0, 0, 0.119748)));
        mol0.addAtom(new Atom("H", new Vector3D(0,  0.761561, -0.478993)));
        mol0.addAtom(new Atom("H", new Vector3D(0, -0.761561, -0.478993)));
        SCFMethod scfm = runSCF(mol0);
        double analyticForce = new HartreeFockForce().computeForce(1, scfm).getY();

        System.out.printf("Water H1 y force: numeric=%.6f  analytic=%.6f%n", numForce, analyticForce);
        assertEquals(numForce, analyticForce, 1e-3);
    }
}
