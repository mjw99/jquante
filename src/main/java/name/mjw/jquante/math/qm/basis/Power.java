package name.mjw.jquante.math.qm.basis;

import net.jafama.FastMath;

/**
 * Represents the powers on orbitals. <br>
 * They are also the magnetic quantum numbers.
 * 
 * @author V.Ganesh
 */
public class Power implements Cloneable {

	/**
	 * Orbital power l.
	 */
	private int l;

	/**
	 * Orbital power m.
	 */
	private int m;

	/**
	 * Orbital power n.
	 */
	private int n;

	/**
	 * Creates a new instance of Power.
	 * 
	 * @param l
	 *            the orbital power l
	 * @param m
	 *            the orbital power m
	 * @param n
	 *            the orbital power n
	 */
	public Power(int l, int m, int n) {
		this.l = l;
		this.m = m;
		this.n = n;
	}

	/**
	 * Getter for property l.
	 * 
	 * @return orbital power l
	 */
	public int getL() {
		return this.l;
	}

	/**
	 * Setter for property l.
	 * 
	 * @param l
	 *            Set orbital power l.
	 */
	public void setL(int l) {
		this.l = l;
	}

	/**
	 * Getter for property m.
	 * 
	 * @return orbital power m.
	 */
	public int getM() {
		return this.m;
	}

	/**
	 * Setter for property m.
	 * 
	 * @param m
	 *            Set orbital power m.
	 */
	public void setM(int m) {
		this.m = m;
	}

	/**
	 * Getter for property n.
	 * 
	 * @return orbital power n.
	 */
	public int getN() {
		return this.n;
	}

	/**
	 * Setter for property n.
	 * 
	 * @param n
	 *            Set orbital power n.
	 */
	public void setN(int n) {
		this.n = n;
	}

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
	 * Clone this !
	 * 
	 * @return the clone
	 */
	@Override
	public Power clone() {
		return new Power(l, m, n);
	}

	/**
	 * overloaded toString()
	 */
	@Override
	public String toString() {
		return "[" + l + ", " + m + ", " + n + "]";
	}
}