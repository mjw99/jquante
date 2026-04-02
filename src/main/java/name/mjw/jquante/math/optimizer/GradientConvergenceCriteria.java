package name.mjw.jquante.math.optimizer;

/**
 * Convergence criterion based on both the energy change and the gradient norm.
 *
 * <p>Convergence is declared when the absolute energy change between steps
 * falls below {@link #tolerance} AND the maximum absolute component of the
 * gradient (or force) falls below {@link #gradientTolerance}.
 */
public class GradientConvergenceCriteria extends ConvergenceCriteria {

    /** Default gradient convergence threshold (Hartree/Bohr). */
    private static final double DEFAULT_GRADIENT_TOLERANCE = 4.5e-4;

    /** The function whose derivatives are checked for convergence. */
    private final OptimizerFunction optimizerFunction;

    /** Maximum absolute gradient component threshold. */
    private double gradientTolerance;

    /**
     * Create a gradient convergence criterion with default tolerances.
     *
     * @param optimizerFunction the function being optimised
     */
    public GradientConvergenceCriteria(OptimizerFunction optimizerFunction) {
        super();
        this.optimizerFunction = optimizerFunction;
        this.gradientTolerance = DEFAULT_GRADIENT_TOLERANCE;
        this.name = "GradientConvergenceCriteria";
    }

    /**
     * Get the gradient convergence tolerance.
     *
     * @return gradient tolerance (Hartree/Bohr)
     */
    public double getGradientTolerance() {
        return gradientTolerance;
    }

    /**
     * Set the gradient convergence tolerance.
     *
     * @param gradientTolerance new tolerance (Hartree/Bohr)
     */
    public void setGradientTolerance(double gradientTolerance) {
        this.gradientTolerance = gradientTolerance;
    }

    /**
     * Return {@code true} when the energy change is below {@link #tolerance}
     * AND the largest gradient component is below {@link #gradientTolerance}.
     *
     * @return {@code true} if converged
     */
    @Override
    public boolean isConverged() {
        boolean energyConverged = Math.abs(newValue - oldValue) < tolerance;
        boolean gradientConverged = optimizerFunction.getMaxNormOfDerivatives() < gradientTolerance;
        return energyConverged && gradientConverged;
    }
}
