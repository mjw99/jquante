package name.mjw.jquante.math.qm.basis;

import net.jafama.FastMath;

/**
 * Represents the powers on orbitals. <br>
 * They are also the magnetic quantum numbers.
 * 
 * @author V.Ganesh
 * @author M.Williamson
 */
public record Power(int l, int m, int n) {

	/**
	 * Return the maximum angular momentum of this Power
	 * 
	 * @return the maximum of the powers
	 */
	public int getMaximumAngularMomentum() {
		return FastMath.max(FastMath.max(l, m), n);
	}

	/**
	 * Return the minimum angular momentum of this Power
	 * 
	 * @return the minimum of the powers
	 */
	public int getMinimumAngularMomentum() {
		return FastMath.min(FastMath.min(l, m), n);
	}

	/**
	 * Return the total angular momentum of this Power
	 * 
	 * @return the maximum of the powers
	 */
	public int getTotalAngularMomentum() {
		return l + m + n;
	}

	/**
	 * Get the power as a dimension 3 array.
	 *
	 * @return powers
	 */
	public double[] toArray() {
		return new double[] { l, m, n };
	}

}
