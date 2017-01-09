package name.mjw.jquante.math.optimizer.impl;

import name.mjw.jquante.math.optimizer.AbstractOptimizer;
import name.mjw.jquante.math.optimizer.OptimizerFunction;

/**
 * The "step it" optimizer.
 * 
 * @author V. Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class StepitOptimizer extends AbstractOptimizer {

	/**
	 * Creates a new instance of StepitOptimizer
	 * 
	 * @param function
	 *            Function to be optmimized.
	 */
	public StepitOptimizer(OptimizerFunction function) {
		super(function);
	}

	/**
	 * Apply the simplex minimizer
	 */
	@Override
	public void minimize() {
		throw new UnsupportedOperationException("Not supported yet.");
	}
}
