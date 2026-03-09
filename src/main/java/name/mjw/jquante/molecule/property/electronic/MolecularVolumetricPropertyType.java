package name.mjw.jquante.molecule.property.electronic;

/**
 * Type of volumetric property.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public enum MolecularVolumetricPropertyType {
	/** Electron density, representing the probability of finding an electron at a given point. */
	ELECTRON_DENSITY,
	/** Molecular electrostatic potential, representing the electrostatic interaction energy with a unit positive charge. */
	ELECTROSTATIC_POTENTIAL,
	/** Electric field, the negative gradient of the electrostatic potential. */
	ELECTRIC_FIELD,
	/** An unrecognised or unspecified volumetric property type. */
	UNKNOWN_TYPE
}
