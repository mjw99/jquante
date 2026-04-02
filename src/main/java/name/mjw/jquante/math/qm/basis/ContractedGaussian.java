package name.mjw.jquante.math.qm.basis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Objects;

import name.mjw.jquante.config.impl.AtomInfo;
import org.hipparchus.geometry.euclidean.threed.Vector3D;

import com.google.common.collect.ComparisonChain;

import name.mjw.jquante.molecule.Atom;
import name.mjw.jquante.molecule.Molecule;
import net.jafama.FastMath;

/**
 * The class defines a contracted Gaussian and the operations on it.
 * 
 * @author V.Ganesh
 * @author mw529
 * @version 2.0 (Part of MeTA v2.0)
 */
public class ContractedGaussian implements Comparable<ContractedGaussian> {

	/**
	 * Holds value of property origin.
	 */
	private final Vector3D origin;

	/**
	 * Holds value of property powers, whole array size is 3.
	 */
	private final Power powers;

	/**
	 * The Primitive Gaussians (PGs) that make up this Contracted Gaussian.
	 */
	private final ArrayList<PrimitiveGaussian> primitives;

	/**
	 * The list of exponents for this Gaussian.
	 */
	private final ArrayList<Double> exponents;

	/**
	 * Holds value of property coefficients.
	 */
	private final ArrayList<Double> coefficients;

	/**
	 * Normalization factors for each PG
	 */
	private ArrayList<Double> primNorms;

	/**
	 * Holds value of property normalization.
	 */
	private double normalization;

	/** The atom on which this contracted Gaussian is centered, or null if centered on a point. */
	protected Atom centeredAtom;

	/** The index of this contracted Gaussian within the overall basis function list. */
	protected int basisFunctionIndex;

	/**
	 * Creates a new instance of ContractedGaussian
	 * 
	 * @param origin - the (x, y, z) on which this gaussian is centered
	 * @param powers - the powers of this gaussian
	 */
	public ContractedGaussian(Vector3D origin, Power powers) {
		this.origin = origin;
		this.powers = powers;

		primitives = new ArrayList<>();
		exponents = new ArrayList<>();
		coefficients = new ArrayList<>();
		primNorms = new ArrayList<>();

		this.normalization = 1;
	}

	/**
	 * Creates a new instance of ContractedGaussian
	 * 
	 * @param centeredAtom - CG centered on an atom
	 * @param powers       - the powers of this gaussian
	 */
	public ContractedGaussian(Atom centeredAtom, Power powers) {
		this(centeredAtom.getAtomCenterInAU(), powers);

		this.centeredAtom = centeredAtom;
	}

	/**
	 * Creates a new instance of ContractedGaussian
	 *
	 * @param cg ContractedGaussian
	 */
	public ContractedGaussian(ContractedGaussian cg) {
		this.origin = cg.getOrigin();
		this.powers = cg.getPowers();

		primitives = cg.getPrimitives();
		exponents = cg.getExponents();
		coefficients = cg.getCoefficients();
		primNorms = cg.getPrimNorms();
	}

	/**
	 * Get the value of centeredAtom
	 * 
	 * @return the value of centeredAtom
	 */
	public Atom getCenteredAtom() {
		return centeredAtom;
	}

	/**
	 * Adds a primitive gaussian (PG) to this contracted gaussian list. This method
	 * will re-sort the primitives in this class according to the
	 * {@link PrimitiveGaussian#compareTo Libint logic}, everytime another primitive
	 * is added.
	 * 
	 * @param exponent    the exponent for this PG
	 * @param coefficient the coefficient of this PG
	 */
	public void addPrimitive(double exponent, double coefficient) {
		primitives.add(new PrimitiveGaussian(origin, powers, exponent, coefficient));

		Collections.sort(primitives);
		primitives.trimToSize();

		this.exponents.clear();
		for (int i = 0; i < this.primitives.size(); i++) {
			exponents.add(primitives.get(i).exponent());
		}
		exponents.trimToSize();

		this.coefficients.clear();
		for (int i = 0; i < this.primitives.size(); i++) {
			coefficients.add(primitives.get(i).coefficient());
		}
		coefficients.trimToSize();

		this.primNorms.clear();
		for (int i = 0; i < this.primitives.size(); i++) {
			primNorms.add(primitives.get(i).normalization());
		}
		primNorms.trimToSize();
	}

	/**
	 * Getter for property origin.
	 * 
	 * @return Value of property origin.
	 */
	public Vector3D getOrigin() {
		return this.origin;
	}

	/**
	 * Getter for property powers.
	 * 
	 * @return Value of property powers.
	 */
	public Power getPowers() {
		return this.powers;
	}

	/**
	 * Getter for property primitives.
	 * 
	 * @return Value of property primitives.
	 */
	public ArrayList<PrimitiveGaussian> getPrimitives() {
		return this.primitives;
	}

	/**
	 * Getter for property exponents.
	 * 
	 * @return Value of property exponents.
	 */
	public ArrayList<Double> getExponents() {
		return this.exponents;
	}

	/**
	 * Getter for property coefficients.
	 * 
	 * @return Value of property coefficients.
	 */
	public ArrayList<Double> getCoefficients() {
		return this.coefficients;
	}

	/**
	 * Getter for property primNorms.
	 *	
	 * @return Value of property primNorms.
	 */
	public ArrayList<Double> getPrimNorms() {
		return this.primNorms;
	}

	/**
	 * Getter for property normalisation.
	 * 
	 * @return Value of property normalisation.
	 */
	public double getNormalization() {
		return this.normalization;
	}

	/**
	 * Normalise this basis function.
	 * 
	 */
	public void normalize() {
		normalization = 1.0 / FastMath.sqrt(this.overlap(this));

		this.primNorms.clear();
		for (int i = 0; i < primitives.size(); i++) {
			primNorms.add(primitives.get(i).normalization());
		}
	}

	/**
	 * Overlap matrix element with another ContractedGaussian
	 * 
	 * @param cg the ContractedGaussian with which the overlap is to be be
	 *           determined.
	 * @return the overlap value
	 */
	public double overlap(ContractedGaussian cg) {
		double sij = 0.0;

		ArrayList<PrimitiveGaussian> cgPrimitives = cg.getPrimitives();

		for (int i = 0; i < primitives.size(); i++) {
			PrimitiveGaussian iPG = primitives.get(i);
			for (int j = 0; j < cgPrimitives.size(); j++) {
				PrimitiveGaussian jPG = cgPrimitives.get(j);

				sij += iPG.coefficient() * jPG.coefficient() * iPG.overlap(jPG);
			}
		}

		return normalization * cg.normalization * sij;
	}

	/**
	 * Compute the overlap energy derivative term w.r.t the specified atom index.
	 *
	 * @param atomIndex the reference atomIndex
	 * @param cg        the other cg
	 * @return the derivative term as an instance of Vector3D
	 */
	public Vector3D overlapDerivative(int atomIndex, ContractedGaussian cg) {
		double dx = 0.0;
		double dy = 0.0;
		double dz = 0.0;

		// case 1: atomIndex is centered on this CG
		if (this.centeredAtom != null && this.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian jPG : cg.primitives) {
				for (PrimitiveGaussian iPG : primitives) {
					Vector3D v = overlapDerivativeHelper(iPG, jPG, origin);
					dx += v.getX(); dy += v.getY(); dz += v.getZ();
				}
			}
		}

		// case 2: atomIndex is centered on the other CG
		if (cg.centeredAtom != null && cg.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian iPG : primitives) {
				for (PrimitiveGaussian jPG : cg.primitives) {
					Vector3D v = overlapDerivativeHelper(jPG, iPG, cg.origin);
					dx += v.getX(); dy += v.getY(); dz += v.getZ();
				}
			}
		}

		// case 3: atomIndex is not centered on either CG → derivative is zero.

		return new Vector3D(normalization * cg.normalization * dx,
				normalization * cg.normalization * dy,
				normalization * cg.normalization * dz);
	}

	/**
	 * Compute the contribution of one primitive pair to the overlap derivative
	 * with respect to the centre {@code pOrigin} of {@code iPG}.
	 *
	 * <p>The formula for the x-component is:
	 * <pre>
	 *   ∂⟨iPG|jPG⟩/∂X = sqrt(α(2l+1))·c·⟨(l+1,m,n)|jPG⟩
	 *                   − 2l·sqrt(α/(2l−1))·c·⟨(l−1,m,n)|jPG⟩
	 * </pre>
	 * and analogously for y (m) and z (n).
	 */
	private Vector3D overlapDerivativeHelper(PrimitiveGaussian iPG,
			PrimitiveGaussian jPG, Vector3D pOrigin) {
		int l = iPG.powers().l();
		int m = iPG.powers().m();
		int n = iPG.powers().n();
		double coeff = iPG.coefficient() * jPG.coefficient();
		double alpha = iPG.exponent();

		// X-component
		double termAx = FastMath.sqrt(alpha * (2.0 * l + 1.0)) * coeff
				* new PrimitiveGaussian(pOrigin, new Power(l + 1, m, n), alpha, 1.0).overlap(jPG);
		double termBx = 0.0;
		if (l > 0) {
			termBx = -2.0 * l * FastMath.sqrt(alpha / (2.0 * l - 1.0)) * coeff
					* new PrimitiveGaussian(pOrigin, new Power(l - 1, m, n), alpha, 1.0).overlap(jPG);
		}

		// Y-component
		double termAy = FastMath.sqrt(alpha * (2.0 * m + 1.0)) * coeff
				* new PrimitiveGaussian(pOrigin, new Power(l, m + 1, n), alpha, 1.0).overlap(jPG);
		double termBy = 0.0;
		if (m > 0) {
			termBy = -2.0 * m * FastMath.sqrt(alpha / (2.0 * m - 1.0)) * coeff
					* new PrimitiveGaussian(pOrigin, new Power(l, m - 1, n), alpha, 1.0).overlap(jPG);
		}

		// Z-component
		double termAz = FastMath.sqrt(alpha * (2.0 * n + 1.0)) * coeff
				* new PrimitiveGaussian(pOrigin, new Power(l, m, n + 1), alpha, 1.0).overlap(jPG);
		double termBz = 0.0;
		if (n > 0) {
			termBz = -2.0 * n * FastMath.sqrt(alpha / (2.0 * n - 1.0)) * coeff
					* new PrimitiveGaussian(pOrigin, new Power(l, m, n - 1), alpha, 1.0).overlap(jPG);
		}

		return new Vector3D(termAx + termBx, termAy + termBy, termAz + termBz);
	}

	/**
	 * Kinetic Energy (KE) matrix element with another ContractedGaussian
	 * 
	 * @param cg the ContractedGaussian with which KE is to be determined.
	 * @return the KE value
	 */
	public double kinetic(ContractedGaussian cg) {
		double tij = 0.0;
		List<PrimitiveGaussian> cgPrimitives = cg.getPrimitives();

		for (int i = 0; i < primitives.size(); i++) {
			PrimitiveGaussian iPG = primitives.get(i);
			for (int j = 0; j < cgPrimitives.size(); j++) {
				PrimitiveGaussian jPG = cgPrimitives.get(j);

				tij += iPG.coefficient() * jPG.coefficient() * iPG.kinetic(jPG);
			}
		}

		return normalization * cg.normalization * tij;
	}

	/**
	 * Compute the kinetic energy derivative term w.r.t the specified atom index
	 * 
	 * @param atomIndex the reference atomIndex
	 * @param cg        the other cg
	 * @return the derivative term as an instance of Vector3D
	 */
	public Vector3D kineticDerivative(int atomIndex, ContractedGaussian cg) {
		double dx = 0.0, dy = 0.0, dz = 0.0;

		// case 1: atomIndex is centered on this CG
		if (this.centeredAtom != null && this.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian jPG : cg.primitives) {
				for (PrimitiveGaussian iPG : primitives) {
					Vector3D v = kineticDerivativeHelper(iPG, jPG, origin);
					dx += v.getX(); dy += v.getY(); dz += v.getZ();
				}
			}
		}

		// case 2: atomIndex is centered on the other CG
		if (cg.centeredAtom != null && cg.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian iPG : primitives) {
				for (PrimitiveGaussian jPG : cg.primitives) {
					Vector3D v = kineticDerivativeHelper(jPG, iPG, cg.origin);
					dx += v.getX(); dy += v.getY(); dz += v.getZ();
				}
			}
		}

		// case 3: atomIndex not centered on either CG → derivative is zero.

		return new Vector3D(normalization * cg.normalization * dx,
				normalization * cg.normalization * dy,
				normalization * cg.normalization * dz);
	}

	/** helper method for kinetic energy derivative */
	private Vector3D kineticDerivativeHelper(PrimitiveGaussian iPG, PrimitiveGaussian jPG, Vector3D pOrigin) {
		int l = iPG.powers().l();
		int m = iPG.powers().m();
		int n = iPG.powers().n();
		double coeff = iPG.coefficient() * jPG.coefficient();
		double alpha = iPG.exponent();

		// X-component
		double termAx = FastMath.sqrt(alpha * (2.0 * l + 1.0)) * coeff
				* new PrimitiveGaussian(pOrigin, new Power(l + 1, m, n), alpha, 1.0).kinetic(jPG);
		double termBx = 0.0;
		if (l > 0) {
			termBx = -2.0 * l * FastMath.sqrt(alpha / (2.0 * l - 1.0)) * coeff
					* new PrimitiveGaussian(pOrigin, new Power(l - 1, m, n), alpha, 1.0).kinetic(jPG);
		}

		// Y-component
		double termAy = FastMath.sqrt(alpha * (2.0 * m + 1.0)) * coeff
				* new PrimitiveGaussian(pOrigin, new Power(l, m + 1, n), alpha, 1.0).kinetic(jPG);
		double termBy = 0.0;
		if (m > 0) {
			termBy = -2.0 * m * FastMath.sqrt(alpha / (2.0 * m - 1.0)) * coeff
					* new PrimitiveGaussian(pOrigin, new Power(l, m - 1, n), alpha, 1.0).kinetic(jPG);
		}

		// Z-component
		double termAz = FastMath.sqrt(alpha * (2.0 * n + 1.0)) * coeff
				* new PrimitiveGaussian(pOrigin, new Power(l, m, n + 1), alpha, 1.0).kinetic(jPG);
		double termBz = 0.0;
		if (n > 0) {
			termBz = -2.0 * n * FastMath.sqrt(alpha / (2.0 * n - 1.0)) * coeff
					* new PrimitiveGaussian(pOrigin, new Power(l, m, n - 1), alpha, 1.0).kinetic(jPG);
		}

		return new Vector3D(termAx + termBx, termAy + termBy, termAz + termBz);
	}

	/**
	 * Nuclear matrix element with another ContractedGaussian and a center
	 * 
	 * @param cg     the ContractedGaussian with which nuclear interaction is to
	 *               determined.
	 * @param center the center at which nuclear energy is to be computed
	 * @return the nuclear value
	 */
	public double nuclear(ContractedGaussian cg, Vector3D center) {
		double vij = 0.0;
		int i;
		int j;

		List<PrimitiveGaussian> cgPrimitives = cg.getPrimitives();

		for (i = 0; i < primitives.size(); i++) {
			PrimitiveGaussian iPG = primitives.get(i);
			for (j = 0; j < cgPrimitives.size(); j++) {
				PrimitiveGaussian jPG = cgPrimitives.get(j);

				vij += iPG.coefficient() * jPG.coefficient() * iPG.nuclear(jPG, center);
			} // end for
		} // end for

		return normalization * cg.normalization * vij;
	}

	/**
	 * Derivative of nuclear attraction term, with respect to an atom index
	 * 
	 * @param mol       the reference Molecule object
	 * @param atomIndex the reference atomIndex
	 * @param cg        the other ContractedGaussian
	 * @return partial derivative of nuclear attraction integral term
	 */
	public Vector3D nuclearAttractionDerivative(Molecule mol, int atomIndex, ContractedGaussian cg) {
		double dx = 0.0, dy = 0.0, dz = 0.0;

		// case 1: atom index is centered on this CG
		if (this.centeredAtom != null && this.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian jPG : cg.primitives) {
				for (PrimitiveGaussian iPG : primitives) {
					Vector3D v = nuclearDerivativeHelper(mol, iPG, jPG, origin);
					dx += v.getX(); dy += v.getY(); dz += v.getZ();
				}
			}
		}

		// case 2: atomIndex is centered on the other CG
		if (cg.centeredAtom != null && cg.centeredAtom.getIndex() == atomIndex) {
			for (PrimitiveGaussian iPG : primitives) {
				for (PrimitiveGaussian jPG : cg.primitives) {
					Vector3D v = nuclearDerivativeHelper(mol, jPG, iPG, cg.origin);
					dx += v.getX(); dy += v.getY(); dz += v.getZ();
				}
			}
		}

		// case 3: derivative w.r.t. the nuclear centre at atomIndex.
		// Z_A * d<μ|1/|r-R_A||ν>/dR_A is nonzero for ALL (μ,ν) pairs regardless
		// of where μ and ν are centred, so this is computed unconditionally.
		{
			Atom atomA = mol.getAtom(atomIndex);
			double atno = AtomInfo.getInstance().getAtomicNumber(atomA.getSymbol());
			for (PrimitiveGaussian iPG : primitives) {
				for (PrimitiveGaussian jPG : cg.primitives) {
					double factor = iPG.coefficient() * jPG.coefficient() * atno;
					Vector3D grad = jPG.nuclearAttractionGradient(iPG, atomA.getAtomCenter()).scalarMultiply(factor);
					dx += grad.getX(); dy += grad.getY(); dz += grad.getZ();
				}
			}
		}

		return new Vector3D(normalization * cg.normalization * dx,
				normalization * cg.normalization * dy,
				normalization * cg.normalization * dz);
	}

	/** helper method for nuclear attraction derivative */
	private Vector3D nuclearDerivativeHelper(Molecule mol, PrimitiveGaussian iPG, PrimitiveGaussian jPG,
			Vector3D pOrigin) {
		int l = iPG.powers().l();
		int m = iPG.powers().m();
		int n = iPG.powers().n();
		double coeff = iPG.coefficient() * jPG.coefficient();
		double alpha = iPG.exponent();

		AtomInfo ai = AtomInfo.getInstance();

		// X-component
		double termAx = 0.0;
		PrimitiveGaussian pgLp1 = new PrimitiveGaussian(pOrigin, new Power(l + 1, m, n), alpha, 1.0);
		for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
			Atom atom = mol.getAtom(i);
			termAx += ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha * (2.0 * l + 1.0)) * coeff
					* pgLp1.nuclear(jPG, atom.getAtomCenter());
		}
		double termBx = 0.0;
		if (l > 0) {
			PrimitiveGaussian pgLm1 = new PrimitiveGaussian(pOrigin, new Power(l - 1, m, n), alpha, 1.0);
			for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
				Atom atom = mol.getAtom(i);
				termBx += -2.0 * l * ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha / (2.0 * l - 1.0))
						* coeff * pgLm1.nuclear(jPG, atom.getAtomCenter());
			}
		}

		// Y-component
		double termAy = 0.0;
		PrimitiveGaussian pgMp1 = new PrimitiveGaussian(pOrigin, new Power(l, m + 1, n), alpha, 1.0);
		for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
			Atom atom = mol.getAtom(i);
			termAy += ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha * (2.0 * m + 1.0)) * coeff
					* pgMp1.nuclear(jPG, atom.getAtomCenter());
		}
		double termBy = 0.0;
		if (m > 0) {
			PrimitiveGaussian pgMm1 = new PrimitiveGaussian(pOrigin, new Power(l, m - 1, n), alpha, 1.0);
			for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
				Atom atom = mol.getAtom(i);
				termBy += -2.0 * m * ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha / (2.0 * m - 1.0))
						* coeff * pgMm1.nuclear(jPG, atom.getAtomCenter());
			}
		}

		// Z-component
		double termAz = 0.0;
		PrimitiveGaussian pgNp1 = new PrimitiveGaussian(pOrigin, new Power(l, m, n + 1), alpha, 1.0);
		for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
			Atom atom = mol.getAtom(i);
			termAz += ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha * (2.0 * n + 1.0)) * coeff
					* pgNp1.nuclear(jPG, atom.getAtomCenter());
		}
		double termBz = 0.0;
		if (n > 0) {
			PrimitiveGaussian pgNm1 = new PrimitiveGaussian(pOrigin, new Power(l, m, n - 1), alpha, 1.0);
			for (int i = 0; i < mol.getNumberOfAtoms(); i++) {
				Atom atom = mol.getAtom(i);
				termBz += -2.0 * n * ai.getAtomicNumber(atom.getSymbol()) * FastMath.sqrt(alpha / (2.0 * n - 1.0))
						* coeff * pgNm1.nuclear(jPG, atom.getAtomCenter());
			}
		}

		return new Vector3D(termAx + termBx, termAy + termBy, termAz + termBz);
	}

	/**
	 * The amplitude of this primitive Gaussian at a given point.
	 * 
	 * @param point the reference point
	 * @return the amplitude of this PG at the specified point
	 */
	public double amplitude(Vector3D point) {
		double value = 0.0;

		for (int i = 0; i < primitives.size(); i++) {
			value += primitives.get(i).amplitude(point);
		}
		return (normalization * value);
	}

	/**
	 * Calculate laplacian at the specified point
	 * 
	 * @param point the point where Laplacian is to be computed
	 * @return the laplacian at this point
	 */
	public double laplacian(Vector3D point) {
		double value = 0.0;

		for (PrimitiveGaussian pg : primitives)
			value += pg.laplacian(point);

		return normalization * value;
	}

	/**
	 * Evaluate gradient of this function at a given point
	 * 
	 * @param point the point where gradient is to be evaluated
	 * @return partial derivatives with respect to x, y, z
	 */
	public Vector3D gradient(Vector3D point) {
		Vector3D grad = new Vector3D(0, 0, 0);

		for (PrimitiveGaussian pg : primitives)
			grad = grad.add(pg.gradient(point));

		return grad.scalarMultiply(normalization);
	}

	/**
	 * Return the maximum angular momentum for this contracted basis function
	 * 
	 * @return the maximum angular momentum of this contracted basis function
	 */
	public int getMaximumAngularMomentum() {
		return powers.getMaximumAngularMomentum();
	}

	/**
	 * Return the minimum angular momentum for this contracted basis function
	 * 
	 * @return the minimum angular momentum of this contracted basis function
	 */
	public int getMinimumAngularMomentum() {
		return powers.getMinimumAngularMomentum();
	}

	/**
	 * Return the total angular momentum of this contracted basis function
	 * 
	 * @return the maximum of the contracted basis function
	 */
	public int getTotalAngularMomentum() {
		return powers.getTotalAngularMomentum();
	}

	/**
	 * Get the value of index
	 * 
	 * @return the value of index
	 */
	public int getBasisFunctionIndex() {
		return basisFunctionIndex;
	}

	/**
	 * Set the value of index
	 * 
	 * @param index new value of index
	 */
	public void setBasisFunctionIndex(int index) {
		this.basisFunctionIndex = index;
	}

	/**
	 * Determines if two ContractedGaussians belong to the same shell.
	 * 
	 * A shell is a set of basis function with the same centre and same contracted
	 * exponents, for example, sp shell, d shell.
	 *
	 * @param that The other contracted Gaussian.
	 * @return True if otherCg is in the same shell.
	 */
	public boolean isSameShell(ContractedGaussian that) {

		if (!this.origin.equals(that.origin)) {
			return false;
		}

		if (this.getPrimitives().size() != that.getPrimitives().size()) {
			return false;
		}

		if (this.getTotalAngularMomentum() != that.getTotalAngularMomentum()) {
			return false;
		}

		for (int i = 0; i < this.getExponents().size(); i++) {
			if (!this.getExponents().get(i).equals(that.getExponents().get(i))) {

				return false;
			}

		}

		for (int i = 0; i < this.getCoefficients().size(); i++) {
			if (!this.getCoefficients().get(i).equals(that.getCoefficients().get(i))) {

				return false;
			}

		}

		return true;
	}

	/**
	 * Returns a flattened 1D double array representation of this contracted Gaussian,
	 * suitable for passing to native/GPU kernels. The layout is:
	 * [x, y, z, l, m, n, nExp, c0..cN, e0..eN, norm0..normN].
	 *
	 * @return a flat double array encoding this contracted Gaussian's parameters
	 */
	public double[] flattern() {

		int preDynamicSize = 3 + 3 + 1;
		int numberOfExp = this.getCoefficients().size();

		double[] flatContractedGaussian = new double[preDynamicSize + numberOfExp + numberOfExp + numberOfExp];

		flatContractedGaussian[0] = this.origin.getX();
		flatContractedGaussian[1] = this.origin.getY();
		flatContractedGaussian[2] = this.origin.getZ();

		flatContractedGaussian[3] = this.powers.l();
		flatContractedGaussian[4] = this.powers.m();
		flatContractedGaussian[5] = this.powers.n();

		flatContractedGaussian[6] = numberOfExp;

		int offset = preDynamicSize;

		for (int i = 0; i < this.getCoefficients().size(); i++) {

			flatContractedGaussian[i + offset] = this.getCoefficients().get(i);
		}

		offset = offset + numberOfExp;
		for (int i = 0; i < this.getCoefficients().size(); i++) {

			flatContractedGaussian[i + offset] = this.getExponents().get(i);
		}

		offset = offset + numberOfExp;
		for (int i = 0; i < this.getCoefficients().size(); i++) {

			flatContractedGaussian[i + offset] = this.getPrimNorms().get(i);
		}

		return flatContractedGaussian;
	}
	/**
	 * overloaded toString()
	 */
	@Override
	public String toString() {
		return "Origin : " + origin + " Powers : " + powers + " " + exponents + " " + coefficients + "\n";
	}

	/**
	 * Returns a hash code for this contracted Gaussian based on all its structural parameters.
	 *
	 * @return a hash code value for this contracted Gaussian
	 */
	@Override
	public int hashCode() {
		return Objects.hash(centeredAtom, coefficients, exponents, basisFunctionIndex, normalization, origin, powers, primNorms,
				primitives);
	}

	/**
	 * Compares this contracted Gaussian to another using the Libint ordering convention:
	 * first by increasing basis function index of the center, then by increasing angular
	 * momentum, and finally by decreasing exponent.
	 *
	 * @param other the other ContractedGaussian to compare to
	 * @return a negative integer, zero, or a positive integer as this object is
	 *         less than, equal to, or greater than the other
	 */
	@Override
	public int compareTo(ContractedGaussian other) {
	    return ComparisonChain.start()
	      .compare( basisFunctionIndex, other.basisFunctionIndex )
	      .compare( powers.getMaximumAngularMomentum(), other.powers.getMaximumAngularMomentum() )
	      .compare( this.powers.l(), other.powers.l() )
	      .compare( this.powers.m(), other.powers.m() )
	      .compare( this.powers.n(), other.powers.n() )
	      .result();
	  }

	/**
	 * Indicates whether some other object is equal to this contracted Gaussian.
	 * Two contracted Gaussians are equal if all their structural parameters match.
	 *
	 * @param obj the reference object with which to compare
	 * @return true if this object is equal to obj, false otherwise
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		ContractedGaussian other = (ContractedGaussian) obj;
		return Objects.equals(centeredAtom, other.centeredAtom) && Objects.equals(coefficients, other.coefficients)
				&& Objects.equals(exponents, other.exponents) && basisFunctionIndex == other.basisFunctionIndex
				&& Double.doubleToLongBits(normalization) == Double.doubleToLongBits(other.normalization)
				&& Objects.equals(origin, other.origin) && Objects.equals(powers, other.powers)
				&& Objects.equals(primNorms, other.primNorms) && Objects.equals(primitives, other.primitives);
	}
}
