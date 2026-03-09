package name.mjw.jquante.common;

/**
 * Enum defining Units used to represent various numbers used in IDE.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public enum Units {
	/** Atomic units (a.u.), the natural unit system for quantum chemistry calculations. */
	AU {
		@Override
		public String toString() {
			return "a.u.";
		}
	},
	/** Angstroms (Å), commonly used for molecular geometry coordinates. */
	ANGSTROM {
		@Override
		public String toString() {
			return "angstroms";
		}
	},
	/** Kilocalories per mole (kcal/mol), used for energy values. */
	KCAL {
		@Override
		public String toString() {
			return "kcal";
		}
	}
}
