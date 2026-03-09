package name.mjw.jquante.molecule;

import java.io.Serializable;
import java.util.EnumSet;

/**
 * This class defines an enum type called BondType, used to differentiate among
 * different classes of chemical bond.
 * 
 * @author V.Ganesh, J. Milthorpe
 * @version 2.0 (Part of MeTA v2.0)
 */
public enum BondType implements Serializable {

	/** A weak, non-covalent interaction such as a van der Waals contact. */
	WEAK_BOND("Weak bond", 0.0),
	/** Absence of any bonding interaction between two atoms. */
	NO_BOND("No bond", 0.0),
	/** A covalent bond involving one shared electron pair. */
	SINGLE_BOND("Single bond", 1.0),
	/** A covalent bond involving two shared electron pairs. */
	DOUBLE_BOND("Double bond", 2.0),
	/** A covalent bond involving three shared electron pairs. */
	TRIPLE_BOND("Triple bond", 3.0),
	/** A covalent bond involving four shared electron pairs. */
	QUADRUPLE_BOND("Quadruple bond", 4.0),
	/** A delocalized bond found in aromatic ring systems, with a formal order of 1.5. */
	AROMATIC_BOND("Aromatic bond", 1.5),
	/** A partially double-bonded resonance structure found in amide linkages, with a formal order of 1.41. */
	AMIDE_BOND("Amide bond", 1.41),
	/** An electrostatic interaction between oppositely charged ions. */
	IONIC_BOND("Ionic bond", 0.0);

	/** Human-readable description of this bond type. */
	private final String description;

	/** The formal bond order associated with this bond type. */
	private final double bondOrder;

	/** Set of all BondType values, used for reverse lookup by description. */
	private static final EnumSet<BondType> vals = EnumSet.allOf(BondType.class);

	/**
	 * Constructs a BondType enum constant with the given description and bond order.
	 *
	 * @param description a human-readable string describing the bond type
	 * @param bondOrder   the formal bond order for this bond type
	 */
	private BondType(String description, double bondOrder) {
		this.description = description;
		this.bondOrder = bondOrder;
	}

	/**
	 * Returns a human-readable description of this bond type.
	 *
	 * @return a string describing this bond type
	 */
	@Override
	public String toString() {
		return description;
	}

	/**
	 * @return the bond order, which in a covalent bond is typically the number
	 *         of valence electrons this atom contributes to the bond.
	 */
	public double getBondOrder() {
		return bondOrder;
	}
	/**
	 * complimentary to toString()
	 * 
	 * @param description
	 *            a string indicating bond type
	 * @return the equivalent BondType object
	 */
	public static BondType getBondTypeFor(String description) {
		for (BondType bondType : vals) {
			if (bondType.description.equals(description)) {
				return bondType;
			}
		}
		return BondType.NO_BOND;
	}

	/**
	 * @return true if this is a strong bond type
	 */
	public boolean isStrongBond() {
		return (this != BondType.NO_BOND && this != BondType.WEAK_BOND);
	}
}
