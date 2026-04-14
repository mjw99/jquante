package name.mjw.jquante.math.interpolater;

/**
 * Similar to {@link TriLinearInterpolater} but uses a configurable sub-interpolater
 * to perform individual point interpolation.
 *
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public class ThreePointInterpolater extends Interpolater {

	/** Creates a new instance of ThreePointInterpolater */
	public ThreePointInterpolater() {
		subInterpolater = new LinearInterpolater();
	}

	/** Interpolated value along the first axis at index 1. */
	private double i1;
	/** Interpolated value along the first axis at index 2. */
	private double i2;
	/** Intermediate interpolated value along the second axis at index 1. */
	private double j1;
	/** Intermediate interpolated value along the second axis at index 2. */
	private double j2;
	/** Weight along the first axis. */
	private double w1;
	/** Weight along the second axis. */
	private double w2;
	/** Grid values along the x-axis (9 points). */
	private double[] x1 = new double[9];
	/** Grid values along the y-axis (8 points). */
	private double[] y1 = new double[8];
	/** Delta offsets in each of the three interpolation dimensions. */
	private double[] delta = new double[3];

	/**
	 * Interpolates a value at the given point using a three-point scheme, delegating
	 * individual axis interpolations to the configured sub-interpolater (default:
	 * {@link LinearInterpolater}).
	 *
	 * <p>Note: the sub-interpolater must be a one-point interpolater; using a
	 * multi-point interpolater will produce incorrect results.
	 *
	 * @param y the Y values (results of function evaluation at grid points)
	 * @param x the X values at which interpolation is performed or expected
	 * @return the interpolated value
	 */
	@Override
	public double interpolate(double[] y, double[] x) {
		delta[0] = (x[3] - x[0]) / x[6];
		delta[1] = (x[4] - x[1]) / x[7];
		delta[2] = (x[5] - x[2]) / x[8];

		y1[0] = y[0];
		y1[1] = y[1];
		y1[2] = y[2];
		y1[3] = y[3];
		x1[0] = delta[2];
		i1 = subInterpolater.interpolate(y1, x1);

		y1[0] = y[2];
		y1[1] = y[3];
		y1[2] = y[4];
		y1[3] = y[5];
		i2 = subInterpolater.interpolate(y1, x1);

		y1[0] = y[4];
		y1[1] = y[5];
		y1[2] = y[6];
		y1[3] = y[7];
		j1 = subInterpolater.interpolate(y1, x1);

		y1[0] = y[6];
		y1[1] = y[7];
		y1[2] = y[4];
		y1[3] = y[5];
		j2 = subInterpolater.interpolate(y1, x1);

		y1[0] = i1;
		y1[1] = i2;
		y1[2] = j1;
		y1[3] = j2;
		x1[0] = delta[1];
		w1 = subInterpolater.interpolate(y1, x1);

		y1[0] = j1;
		y1[1] = j2;
		y1[2] = i1;
		y1[3] = i2;
		w2 = subInterpolater.interpolate(y1, x1);

		y1[0] = w1;
		y1[1] = w2;
		y1[2] = w1;
		y1[3] = w2;
		x1[0] = delta[0];

		return (subInterpolater.interpolate(y1, x1));
	}

	/**
	 * Return number of Y arguments required for this Interpolater
	 * 
	 * @return number of Y arguments
	 */
	@Override
	public int getNumberOfRequiredYArgs() {
		return 8;
	}

	/**
	 * Return number of X arguments required for this Interpolater
	 * 
	 * @return number of X arguments
	 */
	@Override
	public int getNumberOfRequiredXArgs() {
		return 9;
	}
}