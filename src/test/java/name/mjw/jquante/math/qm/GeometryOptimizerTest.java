package name.mjw.jquante.math.qm;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.common.Utility;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import name.mjw.jquante.molecule.impl.MoleculeImpl;

class GeometryOptimizerTest {

    private static final Logger LOG = LogManager.getLogger(GeometryOptimizerTest.class);

    /**
     * Optimise H₂ with STO-3G starting from a stretched geometry (0.85 Å)
     * and verify convergence toward the known equilibrium bond length (~0.733 Å).
     *
     * <p>The HF/STO-3G equilibrium bond length of H₂ is 1.3851 bohr = 0.7330 Å.
     */
    @Test
    void optimiseHydrogenSTO3G() throws Exception {
        // Start from a stretched H₂ geometry (longer than equilibrium 0.74 Å).
        Atom h1 = new Atom("H", new Vector3D(0.0, 0.0, 0.0));
        Atom h2 = new Atom("H", new Vector3D(0.85, 0.0, 0.0));

        Molecule hydrogen = new MoleculeImpl("hydrogen");
        hydrogen.addAtom(h1);
        hydrogen.addAtom(h2);

        GeometryOptimizer opt = new GeometryOptimizer(hydrogen, "sto-3g", SCFType.HARTREE_FOCK);
        opt.setInitialStepSize(0.1);
        opt.setGradientThreshold(1e-4);
        opt.setEnergyThreshold(1e-6);
        opt.setMaxIterations(100);

        opt.optimize();

        LOG.info("Converged: {}  iterations: {}  energy: {}",
                opt.isConverged(), opt.getIterations(), opt.getOptimizedEnergy());

        // Record final bond length in Bohr.
        Vector3D r1 = hydrogen.getAtom(0).getAtomCenterInAU();
        Vector3D r2 = hydrogen.getAtom(1).getAtomCenterInAU();
        double bondLengthBohr = r1.distance(r2);
        LOG.info("Final bond length: {} bohr ({} Å)",
                bondLengthBohr, bondLengthBohr * Utility.AU_TO_ANGSTROM_FACTOR);

        assertTrue(opt.isConverged(), "Geometry optimisation should converge.");

        // Energy should be lower than the single-point at 0.85 Å.
        double singlePointEnergy = singlePointEnergy(0.85);
        assertTrue(opt.getOptimizedEnergy() < singlePointEnergy,
                "Optimised energy should be lower than at starting geometry.");

        // Bond length should be near the HF/STO-3G equilibrium value of ~1.385 bohr.
        assertEquals(1.385, bondLengthBohr, 0.05,
                "Optimised bond length should be close to the HF/STO-3G equilibrium (~1.385 bohr).");
    }

    /**
     * Verify that the numerical gradient on H₂ at a stretched geometry points
     * toward contraction — i.e. the gradient component along the bond axis for
     * H1 is negative (moving H1 in -x shortens the bond) and for H2 is positive
     * (moving H2 in +x shortens the bond). Under steepest descent
     * {@code R_new = R_old - step * gradient}, these signs drive atoms together.
     */
    @Test
    void gradientsPointTowardEquilibrium() throws Exception {
        // H₂ at 0.85 Å — longer than equilibrium (~0.74 Å), so gradient should
        // shorten the bond.
        Atom h1 = new Atom("H", new Vector3D(0.0, 0.0, 0.0));
        Atom h2 = new Atom("H", new Vector3D(0.85, 0.0, 0.0));

        Molecule hydrogen = new MoleculeImpl("hydrogen");
        hydrogen.addAtom(h1);
        hydrogen.addAtom(h2);

        GeometryOptimizer opt = new GeometryOptimizer(hydrogen, "sto-3g", SCFType.HARTREE_FOCK);
        double[] gradient = opt.numericalGradient();

        // gradient is [gx_H1, gy_H1, gz_H1, gx_H2, gy_H2, gz_H2] in Hartree/Bohr.
        // At stretched geometry the energy decreases as atoms approach, so:
        //   gx_H1 < 0  (move H1 toward +x to lower energy → -step*grad moves it right)
        //   gx_H2 > 0  (move H2 toward -x to lower energy → -step*grad moves it left)
        double gx_h1 = gradient[0];
        double gx_h2 = gradient[3];

        LOG.info("Numerical gradient: H1=({},{},{})  H2=({},{},{})",
                gradient[0], gradient[1], gradient[2],
                gradient[3], gradient[4], gradient[5]);

        assertTrue(gx_h1 < 0.0,
                "Gradient on H1 should be negative along bond axis (H1 moves toward H2).");
        assertTrue(gx_h2 > 0.0,
                "Gradient on H2 should be positive along bond axis (H2 moves toward H1).");
    }

    /** Compute the single-point HF/STO-3G energy for H₂ at the given bond length (Å). */
    private double singlePointEnergy(double bondLengthAngstrom) throws Exception {
        Atom h1 = new Atom("H", new Vector3D(0.0, 0.0, 0.0));
        Atom h2 = new Atom("H", new Vector3D(bondLengthAngstrom, 0.0, 0.0));

        Molecule hydrogen = new MoleculeImpl("hydrogen");
        hydrogen.addAtom(h1);
        hydrogen.addAtom(h2);

        BasisSetLibrary bsl = new BasisSetLibrary(hydrogen, "sto-3g");
        OneElectronIntegrals e1 = new OneElectronIntegrals(bsl, hydrogen);
        TwoElectronIntegrals e2 = new TwoElectronIntegrals(bsl);

        SCFMethod scf = SCFMethodFactory.getInstance()
                .getSCFMethod(hydrogen, e1, e2, SCFType.HARTREE_FOCK);
        scf.scf();
        return scf.getEnergy();
    }
}
