package name.mjw.jquante.math.qm.integral;

/**
 * Enumeration of Integral package types available in MeTA Studio.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public enum IntegralPackageType {
	/** Nuclear attraction integral term. */
	NUCLEAR_TERM,
	/** One-electron integral terms (kinetic energy and overlap). */
	ONE_ELECTRON_TERM,
	/** Two-electron repulsion integrals evaluated using the Huzinaga scheme. */
	TWO_ELECTRON_HUZINAGA,
	/** Two-electron repulsion integrals evaluated using Rys polynomials. */
	TWO_ELECTRON_RYS,
	/** Two-electron repulsion integrals evaluated using the Head-Gordon/Pople scheme. */
	TWO_ELECTRON_HGP
}
