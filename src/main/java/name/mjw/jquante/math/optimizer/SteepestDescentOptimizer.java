package name.mjw.jquante.math.optimizer;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;

/**
 * A simple steepest-descent (gradient-descent) minimiser with a fixed step
 * size.
 *
 * <p>At each iteration the variables are updated as:
 * <pre>
 *   x_new = x_old + stepSize * derivatives
 * </pre>
 * where {@code derivatives} are the values returned by
 * {@link OptimizerFunction#getDerivatives()} (conventionally forces,
 * i.e. −∂E/∂x).
 *
 * <p>Convergence is checked via the supplied {@link ConvergenceCriteria}.
 * If none is set the optimizer runs for the full {@link #maxIterations}.
 */
public class SteepestDescentOptimizer extends AbstractOptimizer {

    private static final Logger LOG = LogManager.getLogger(SteepestDescentOptimizer.class);

    /** Step size applied to the derivative vector each iteration. */
    private double stepSize;

    /**
     * Create a steepest-descent optimiser.
     *
     * @param optimizerFunction the function to minimise
     * @param stepSize          the step size applied to the derivative each iteration
     */
    public SteepestDescentOptimizer(OptimizerFunction optimizerFunction, double stepSize) {
        super(optimizerFunction);
        this.stepSize = stepSize;
    }

    /**
     * Get the step size.
     *
     * @return the current step size
     */
    public double getStepSize() {
        return stepSize;
    }

    /**
     * Set the step size.
     *
     * @param stepSize new step size
     */
    public void setStepSize(double stepSize) {
        this.stepSize = stepSize;
    }

    /**
     * Run the steepest-descent minimisation.
     *
     * <p>The caller must set {@link #variables} before calling this method.
     */
    @Override
    public void minimize() {
        if (variables == null) {
            throw new IllegalStateException("Variables must be set before calling minimize().");
        }

        currentMinima = optimizerFunction.evaluate(variables);
        if (convergenceCriteria != null) {
            convergenceCriteria.setOldValue(currentMinima);
        }

        for (int iteration = 0; iteration < maxIterations; iteration++) {
            double[] derivs = optimizerFunction.getDerivatives();

            // x_new = x_old + stepSize * derivs  (derivs = forces = −gradient)
            for (int i = 0; i < variables.length; i++) {
                variables[i] += stepSize * derivs[i];
            }

            currentMinima = optimizerFunction.evaluate(variables);
            LOG.debug("Iteration {}: minima = {}", iteration, currentMinima);

            if (convergenceCriteria != null) {
                convergenceCriteria.setNewValue(currentMinima);
                if (convergenceCriteria.isConverged()) {
                    LOG.info("Converged after {} iterations.", iteration + 1);
                    break;
                }
                convergenceCriteria.setOldValue(currentMinima);
            }
        }
    }
}
