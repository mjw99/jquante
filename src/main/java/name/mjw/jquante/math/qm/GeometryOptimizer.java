package name.mjw.jquante.math.qm;

import java.util.ArrayList;
import java.util.List;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.hipparchus.geometry.euclidean.threed.Vector3D;

import name.mjw.jquante.common.Utility;
import name.mjw.jquante.math.qm.basis.BasisSetLibrary;
import name.mjw.jquante.molecule.Molecule;

/**
 * Geometry optimiser for the Hartree-Fock method using steepest descent with a
 * backtracking line search.
 *
 * <p>Energy gradients are computed via central finite differences — each
 * Cartesian coordinate is displaced by ±{@link #displacementBohr} to obtain
 * ∂E/∂X_i ≈ (E(+h) − E(−h)) / (2h).  Although more expensive than analytic
 * gradients (2×3N SCF evaluations per geometry step), this is exact to
 * numerical precision and does not depend on derivative integral implementations.
 *
 * <p>At each step the basis set and electron integrals are rebuilt for the
 * current geometry so that the basis function centres are always consistent
 * with the nuclear positions.
 *
 * <p>Convergence is declared when both the largest Cartesian gradient
 * component and the energy change fall below their respective thresholds.
 */
public class GeometryOptimizer {

    private static final Logger LOG = LogManager.getLogger(GeometryOptimizer.class);

    private final Molecule molecule;
    private final String basisName;
    private final SCFType scfType;

    /** Displacement in Bohr used for central-difference gradient evaluation. */
    private double displacementBohr = 0.001;
    /** Minimum step size in Bohr; backtracking halts and aborts the step when reached. */
    private static final double MIN_STEP_SIZE = 1e-10;
    /** Initial step size in Bohr applied to each gradient component. */
    private double initialStepSize = 0.1;
    /** Convergence threshold for the largest Cartesian gradient (Hartree/Bohr). */
    private double gradientThreshold = 1e-4;
    /** Convergence threshold for the energy change between steps (Hartree). */
    private double energyThreshold = 1e-6;
    /** Maximum number of geometry steps. */
    private int maxIterations = 50;

    private double optimizedEnergy = Double.NaN;
    private boolean converged = false;
    private int iterations = 0;
    private final List<Double> energyHistory = new ArrayList<>();

    /**
     * Create a geometry optimiser.
     *
     * @param molecule  the molecule whose geometry will be optimised in-place
     * @param basisName the basis set name (e.g. {@code "sto-3g"})
     * @param scfType   the SCF method to use
     */
    public GeometryOptimizer(Molecule molecule, String basisName, SCFType scfType) {
        this.molecule = molecule;
        this.basisName = basisName;
        this.scfType = scfType;
    }

    /**
     * Run the geometry optimisation. On return the molecule's atom positions
     * will reflect the optimised (or last-step) geometry.
     *
     * @throws Exception if the basis set cannot be loaded
     */
    public void optimize() throws Exception {
        converged = false;
        energyHistory.clear();

        double oldEnergy = Double.MAX_VALUE;
        double stepSize = initialStepSize;
        boolean justBacktracked = false;

        for (iterations = 0; iterations < maxIterations; iterations++) {
            // Energy at the current geometry.
            double energy = scfEnergy(molecule);
            optimizedEnergy = energy;
            energyHistory.add(energy);

            LOG.info("Geometry step {}: E = {}  stepSize = {}", iterations, energy, stepSize);

            // Backtracking: if energy increased after the previous step, halve
            // the step size and restore the previous geometry.
            if (iterations > 0 && energy > oldEnergy) {
                stepSize *= 0.5;
                LOG.debug("Energy increased; step size reduced to {}", stepSize);
                if (stepSize < MIN_STEP_SIZE) {
                    LOG.warn("Step size fell below minimum ({}); stopping backtracking.", MIN_STEP_SIZE);
                    break;
                }
                molecule.resetAtomCoordinates(previousCoordsAng, false);
                iterations--;
                energyHistory.remove(energyHistory.size() - 1);
                oldEnergy = energyHistory.isEmpty() ? Double.MAX_VALUE
                        : energyHistory.get(energyHistory.size() - 1);
                optimizedEnergy = oldEnergy;
                justBacktracked = true;
                continue;
            }

            // Numerical central-difference gradient.
            double[] gradient = numericalGradient();

            double maxGrad = maxAbsComponent(gradient);
            double deltaE = Math.abs(energy - oldEnergy);
            LOG.info("  maxGrad = {}  deltaE = {}", maxGrad, deltaE);

            // After a backtrack the energy is re-evaluated at the restored geometry,
            // making deltaE trivially ~0.  Skip the energy criterion that iteration
            // so only the gradient drives convergence.
            if (!justBacktracked && maxGrad < gradientThreshold && deltaE < energyThreshold) {
                converged = true;
                break;
            }
            justBacktracked = false;

            oldEnergy = energy;
            previousCoordsAng = currentCoordsAng();

            // Steepest-descent step in AU: R_new = R_old − stepSize * gradient.
            int nAtoms = molecule.getNumberOfAtoms();
            double[] newCoordsAng = new double[nAtoms * 3];
            int idx = 0;
            for (int i = 0; i < nAtoms; i++) {
                Vector3D posAU = molecule.getAtom(i).getAtomCenterInAU();
                newCoordsAng[idx]     = (posAU.getX() - stepSize * gradient[idx])
                        * Utility.AU_TO_ANGSTROM_FACTOR;
                newCoordsAng[idx + 1] = (posAU.getY() - stepSize * gradient[idx + 1])
                        * Utility.AU_TO_ANGSTROM_FACTOR;
                newCoordsAng[idx + 2] = (posAU.getZ() - stepSize * gradient[idx + 2])
                        * Utility.AU_TO_ANGSTROM_FACTOR;
                idx += 3;
            }
            molecule.resetAtomCoordinates(newCoordsAng, false);
        }

        if (!converged) {
            LOG.warn("Geometry optimisation did not converge in {} steps.", maxIterations);
        }
    }

    /**
     * Compute the energy gradient via central finite differences.
     *
     * <p>Returns ∂E/∂X_i in Hartree/Bohr for each Cartesian degree of
     * freedom, packed as [∂E/∂X0, ∂E/∂Y0, ∂E/∂Z0, ∂E/∂X1, ...].
     *
     * @return gradient vector (Hartree/Bohr)
     * @throws Exception if the basis set cannot be loaded
     */
    public double[] numericalGradient() throws Exception {
        int nAtoms = molecule.getNumberOfAtoms();
        double[] gradient = new double[nAtoms * 3];
        double h = displacementBohr;

        // Save the current (unperturbed) Angstrom coordinates.
        double[] origCoordsAng = currentCoordsAng();

        int idx = 0;
        for (int atom = 0; atom < nAtoms; atom++) {
            for (int coord = 0; coord < 3; coord++) {
                // +h displacement in AU → Angstrom
                double[] plus = origCoordsAng.clone();
                plus[idx] += h * Utility.AU_TO_ANGSTROM_FACTOR;
                molecule.resetAtomCoordinates(plus, false);
                double ePlus = scfEnergy(molecule);

                // −h displacement
                double[] minus = origCoordsAng.clone();
                minus[idx] -= h * Utility.AU_TO_ANGSTROM_FACTOR;
                molecule.resetAtomCoordinates(minus, false);
                double eMinus = scfEnergy(molecule);

                gradient[idx] = (ePlus - eMinus) / (2.0 * h);
                idx++;
            }
        }

        // Restore the original geometry.
        molecule.resetAtomCoordinates(origCoordsAng, false);
        return gradient;
    }

    // -------------------------------------------------------------------------
    // Helpers
    // -------------------------------------------------------------------------

    /** Run an SCF calculation for the molecule at its current geometry. */
    private double scfEnergy(Molecule mol) throws Exception {
        BasisSetLibrary bsl = new BasisSetLibrary(mol, basisName);
        OneElectronIntegrals oneEI = new OneElectronIntegrals(bsl, mol);
        TwoElectronIntegrals twoEI = new TwoElectronIntegrals(bsl);
        SCFMethod scf = SCFMethodFactory.getInstance()
                .getSCFMethod(mol, oneEI, twoEI, scfType);
        scf.scf();
        return scf.getEnergy();
    }

    /** Snapshot coordinates in Angstrom before the last step. */
    private double[] previousCoordsAng;

    private double[] currentCoordsAng() {
        int nAtoms = molecule.getNumberOfAtoms();
        double[] coords = new double[nAtoms * 3];
        int idx = 0;
        for (int i = 0; i < nAtoms; i++) {
            Vector3D pos = molecule.getAtom(i).getAtomCenter();
            coords[idx]     = pos.getX();
            coords[idx + 1] = pos.getY();
            coords[idx + 2] = pos.getZ();
            idx += 3;
        }
        return coords;
    }

    private static double maxAbsComponent(double[] v) {
        double max = 0.0;
        for (double d : v) {
            double abs = Math.abs(d);
            if (abs > max) max = abs;
        }
        return max;
    }

    // -------------------------------------------------------------------------
    // Getters / setters
    // -------------------------------------------------------------------------

    /** @return the HF energy at the last (or converged) geometry */
    public double getOptimizedEnergy() {
        return optimizedEnergy;
    }

    /** @return {@code true} if the optimisation converged */
    public boolean isConverged() {
        return converged;
    }

    /** @return the number of geometry steps taken */
    public int getIterations() {
        return iterations;
    }

    /** @return energy at each geometry step */
    public List<Double> getEnergyHistory() {
        return energyHistory;
    }

    /**
     * Set the finite-difference displacement (Bohr; default 0.001).
     *
     * @param displacementBohr displacement for central differences
     */
    public void setDisplacementBohr(double displacementBohr) {
        this.displacementBohr = displacementBohr;
    }

    /**
     * Set the initial steepest-descent step size (Bohr; default 0.1).
     *
     * @param stepSize step size in Bohr
     */
    public void setInitialStepSize(double stepSize) {
        this.initialStepSize = stepSize;
    }

    /**
     * Set the gradient convergence threshold (Hartree/Bohr; default 1e-4).
     *
     * @param threshold gradient convergence threshold
     */
    public void setGradientThreshold(double threshold) {
        this.gradientThreshold = threshold;
    }

    /**
     * Set the energy convergence threshold (Hartree; default 1e-6).
     *
     * @param threshold energy convergence threshold
     */
    public void setEnergyThreshold(double threshold) {
        this.energyThreshold = threshold;
    }

    /**
     * Set the maximum number of geometry steps (default 50).
     *
     * @param maxIterations maximum steps
     */
    public void setMaxIterations(int maxIterations) {
        this.maxIterations = maxIterations;
    }
}
