package name.mjw.jquante;

public class Constants {

	// Misc units
	// ===========

	// Speed of light in m/s
	public static final double CLIGHT = 2.99792458e8;
	// Boltzmann constant
	public static final double KBOLTZ = 3.166830e-6;
	// Coulomb's law coeff if R in \AA and resulting E in eV
	public static final double E2 = 14.399;
	// Planck's constant, in Js
	public static final double PLANCK = 6.6260755e-34;

	// Distance units
	// ===============

	// Conversion of length from bohr to angstrom
	public static final double BOHR2ANG = 0.529177249;
	public static final double ANG2BOHR = 1 / BOHR2ANG;

	// Energy units
	// ============

	// Hartree to kcal/mol conversion
	public static final double HARTREE2KCAL = 627.5095;
	public static final double KCAL2HARTREE = 1 / HARTREE2KCAL;

	// Conversion of energy in eV to energy in kcal/mol
	public static final double EV2KCAL = 23.061;
	public static final double KCAL2EV = 1 / EV2KCAL;

	// Hartree to Joule conversion factor
	public static final double HARTREE2JOULE = 4.3597482e-18;
	public static final double JOULE2HARTREE = 1 / HARTREE2JOULE;

	public static final double EV2HARTREE = HARTREE2KCAL / EV2KCAL;
	public static final double HARTREE2EV = 1 / EV2HARTREE;

	// Mass units
	// ===========

	// Conversion from mass in amu to mass in au (m_e)
	public static final double AMU2ME = 1822.882;
	// # Conversion from mass in au (m_e) to mass in amu
	public static final double ME2AMU = 1 / AMU2ME;

	// Time units
	// ===========

	// Conversion from time in au to time in ps
	public static final double TAU2PS = 41341.447;
	public static final double PS2TAU = 1 / TAU2PS;

	// Derived quantities

	// gas constant R = 1.98722 cal/mole/K
	public static final double RGAS = KBOLTZ * HARTREE2KCAL * 1000.0;

}
