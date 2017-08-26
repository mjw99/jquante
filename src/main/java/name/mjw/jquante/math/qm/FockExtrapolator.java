package name.mjw.jquante.math.qm;

/**
 * Generic interface for Fock matrix extrapolation
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public interface FockExtrapolator {

	/** Initialise this interpolator */
	public void init();

	/**
	 * Get the next extrapolated fock matrix
	 * 
	 * @param currentFock
	 *            the current fock
	 * @param overlap
	 *            the overlap matrix
	 * @param density
	 *            the current density matrix
	 * @return extrapolated fock
	 */
	public Fock next(Fock currentFock, Overlap overlap, Density density);
}
