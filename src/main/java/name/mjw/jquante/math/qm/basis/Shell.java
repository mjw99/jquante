package name.mjw.jquante.math.qm.basis;

/**
 * A shell is a set of basis function with the same centre and same contracted
 * exponent, for example, sp shell, d shell.
 * <p>
 * It also contains a mapping back to the original BasisFunctions that it is
 * composed of.
 * 
 */
public final class Shell extends ContractedGaussian {

	/**
	 * Creates a new Shell from an existing ContractedGaussian.
	 *
	 * @param cg the ContractedGaussian that defines this shell's parameters
	 */
	public Shell(ContractedGaussian cg) {
		super(cg);
	}

	/** The index of the first basis function belonging to this shell. */
	private int firstBasisFunctionIndex;

	/** The index of the last basis function belonging to this shell. */
	private int lastBasisFunctionIndex;

	/**
	 * Returns the index of the first basis function in this shell.
	 *
	 * @return the first basis function index
	 */
	public int getFirstBasisFunctionIndex() {
		return firstBasisFunctionIndex;
	}

	/**
	 * Sets the index of the first basis function in this shell.
	 *
	 * @param firstBasisFunctionIndex the first basis function index
	 */
	public void setFirstBasisFunctionIndex(int firstBasisFunctionIndex) {
		this.firstBasisFunctionIndex = firstBasisFunctionIndex;
	}

	/**
	 * Returns the index of the last basis function in this shell.
	 *
	 * @return the last basis function index
	 */
	public int getLastBasisFunctionIndex() {
		return lastBasisFunctionIndex;
	}

	/**
	 * Sets the index of the last basis function in this shell.
	 *
	 * @param lastBasisFunctionIndex the last basis function index
	 */
	public void setLastBasisFunctionIndex(int lastBasisFunctionIndex) {
		this.lastBasisFunctionIndex = lastBasisFunctionIndex;
	}

	/**
	 * Returns a string representation of this shell as a range of basis function indices.
	 *
	 * @return a string of the form "first-last"
	 */
	@Override
	public String toString() {
		return Integer.toString(firstBasisFunctionIndex) + "-" + Integer.toString(lastBasisFunctionIndex);
	}

	/**
	 * Indicates whether some other object is equal to this Shell.
	 * Two shells are equal if they represent the same shell (same center, exponents, and angular momentum).
	 *
	 * @param that the reference object with which to compare
	 * @return true if the two shells are equivalent
	 */
	@Override
	public boolean equals(Object that) {
		if (this == that)
			return true;
		if (that == null)
			return false;
		if (getClass() != that.getClass())
			return false;
		Shell other = (Shell) that;
		return this.isSameShell(other);
	}

}
