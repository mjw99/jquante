package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.test.Fixtures;

/**
 * Tests for HartreeFockForce (analytic HF gradient).
 *
 * Reference values computed with NWChem 7.2.3, RHF.
 * NWChem reports dE/dR (gradient); jquante returns force = -dE/dR.
 *
 * NWChem reorients water so that the molecular plane contains x and z.
 * In jquante the water geometry uses y and z, so: NWChem-x ↔ -jquante-y.
 */
class HartreeFockForceTest {

    // NWChem gradient precision is ~1e-6 Hartree/Bohr; use 1e-4 tolerance
    private final double delta = 1e-4;

    private SCFMethod runSCF(Molecule mol, String basis) {
        BasisSetLibrary bsl = null;
        try {
            bsl = new BasisSetLibrary(mol, basis);
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
    // H2 / STO-3G
    // NWChem gradients (Hartree/Bohr):
    //   H1: dE/dx = -0.027680   force_x = +0.027680
    //   H2: dE/dx = +0.027680   force_x = -0.027680
    // -----------------------------------------------------------------------

    @Test
    void hydrogenH1ForceSTO3G() {
        SCFMethod scfm = runSCF(Fixtures.getHydrogen(), "sto-3g");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(0, scfm);

        assertEquals( 0.027680, f.getX(), delta);
        assertEquals( 0.0,      f.getY(), delta);
        assertEquals( 0.0,      f.getZ(), delta);
    }

    @Test
    void hydrogenH2ForceSTO3G() {
        SCFMethod scfm = runSCF(Fixtures.getHydrogen(), "sto-3g");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(1, scfm);

        assertEquals(-0.027680, f.getX(), delta);
        assertEquals( 0.0,      f.getY(), delta);
        assertEquals( 0.0,      f.getZ(), delta);
    }

    // -----------------------------------------------------------------------
    // HF / STO-3G
    // NWChem gradients (Hartree/Bohr):
    //   H: dE/dx = +0.059082   force_x = -0.059082
    //   F: dE/dx = -0.059082   force_x = +0.059082
    // -----------------------------------------------------------------------

    @Test
    void hydrogenFluorideHForceSTO3G() {
        SCFMethod scfm = runSCF(Fixtures.getHydrogenFluoride(), "sto-3g");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(0, scfm);

        assertEquals(-0.059082, f.getX(), delta);
        assertEquals( 0.0,      f.getY(), delta);
        assertEquals( 0.0,      f.getZ(), delta);
    }

    @Test
    void hydrogenFluorideFForceSTO3G() {
        SCFMethod scfm = runSCF(Fixtures.getHydrogenFluoride(), "sto-3g");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(1, scfm);

        assertEquals( 0.059082, f.getX(), delta);
        assertEquals( 0.0,      f.getY(), delta);
        assertEquals( 0.0,      f.getZ(), delta);
    }

    // -----------------------------------------------------------------------
    // Water / STO-3G
    // NWChem gradients (Hartree/Bohr):
    //   O:  dE/dz = -0.041936   force_z = +0.041936
    //   H1: dE/dy = -0.013109,  dE/dz = +0.020968
    //       force_y = +0.013109, force_z = -0.020968
    //   H2: dE/dy = +0.013109,  dE/dz = +0.020968
    //       force_y = -0.013109, force_z = -0.020968
    // -----------------------------------------------------------------------

    @Test
    void waterOForceSTO3G() {
        SCFMethod scfm = runSCF(Fixtures.getWater(), "sto-3g");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(0, scfm);

        assertEquals( 0.0,      f.getX(), delta);
        assertEquals( 0.0,      f.getY(), delta);
        assertEquals( 0.041936, f.getZ(), delta);
    }

    @Test
    void waterH1ForceSTO3G() {
        SCFMethod scfm = runSCF(Fixtures.getWater(), "sto-3g");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(1, scfm);

        assertEquals( 0.0,      f.getX(), delta);
        assertEquals( 0.013109, f.getY(), delta);
        assertEquals(-0.020968, f.getZ(), delta);
    }

    @Test
    void waterH2ForceSTO3G() {
        SCFMethod scfm = runSCF(Fixtures.getWater(), "sto-3g");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(2, scfm);

        assertEquals( 0.0,      f.getX(), delta);
        assertEquals(-0.013109, f.getY(), delta);
        assertEquals(-0.020968, f.getZ(), delta);
    }


    // -----------------------------------------------------------------------
    // H2 / 6-31G**
    // NWChem gradients (Hartree/Bohr):
    //   H1: dE/dz = -0.005618   force_x = +0.005618
    //   H2: dE/dz = +0.005618   force_x = -0.005618
    // -----------------------------------------------------------------------

    @Test
    void hydrogenH1Force631Gss() {
        SCFMethod scfm = runSCF(Fixtures.getHydrogen(), "6-31gss");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(0, scfm);

        assertEquals( 0.005618, f.getX(), delta);
        assertEquals( 0.0,      f.getY(), delta);
        assertEquals( 0.0,      f.getZ(), delta);
    }

    @Test
    void hydrogenH2Force631Gss() {
        SCFMethod scfm = runSCF(Fixtures.getHydrogen(), "6-31gss");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(1, scfm);

        assertEquals(-0.005618, f.getX(), delta);
        assertEquals( 0.0,      f.getY(), delta);
        assertEquals( 0.0,      f.getZ(), delta);
    }

   

    // -----------------------------------------------------------------------
    // Water / 6-31G**
    // NWChem gradients (Hartree/Bohr; NWChem x ↔ -jquante y):
    //   O:  dE/dz = +0.036847   force_z = -0.036847
    //   H1: dE/dx = -0.019583,  dE/dz = -0.018423
    //       force_y = -0.019583, force_z = +0.018423
    //   H2: force_y = +0.019583, force_z = +0.018423
    // -----------------------------------------------------------------------

    @Test
    void waterOForce631Gss() {
        SCFMethod scfm = runSCF(Fixtures.getWater(), "6-31gss");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(0, scfm);

        assertEquals( 0.0,      f.getX(), delta);
        assertEquals( 0.0,      f.getY(), delta);
        assertEquals(-0.036847, f.getZ(), delta);
    }

    @Test
    void waterH1Force631Gss() {
        SCFMethod scfm = runSCF(Fixtures.getWater(), "6-31gss");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(1, scfm);

        assertEquals( 0.0,      f.getX(), delta);
        assertEquals(-0.019583, f.getY(), delta);
        assertEquals( 0.018423, f.getZ(), delta);
    }

    @Test
    void waterH2Force631Gss() {
        SCFMethod scfm = runSCF(Fixtures.getWater(), "6-31gss");
        HartreeFockForce hff = new HartreeFockForce();
        Vector3D f = hff.computeForce(2, scfm);

        assertEquals( 0.0,      f.getX(), delta);
        assertEquals( 0.019583, f.getY(), delta);
        assertEquals( 0.018423, f.getZ(), delta);
    }
}
