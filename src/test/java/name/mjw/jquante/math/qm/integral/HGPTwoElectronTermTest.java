package name.mjw.jquante.math.qm.integral;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;

import name.mjw.jquante.math.qm.basis.ContractedGaussian;
import name.mjw.jquante.math.qm.basis.Power;
import name.mjw.jquante.test.Fixtures;

class HGPTwoElectronTermTest {
	private final double delta = 0.000001;

	HGPTwoElectronTerm e2 = new HGPTwoElectronTerm();

	static ContractedGaussian cgtoS0 = Fixtures.getCgtoS0();
	static ContractedGaussian cgtoS1 = Fixtures.getCgtoS1();
	static ContractedGaussian cgtoP0 = Fixtures.getCgtoP0();
	static ContractedGaussian cgtoD0 = Fixtures.getCgtoD0();
	static ContractedGaussian cgtoF0 = Fixtures.getCgtoF0();

	@Test
	void ss00() {
		assertEquals(1.1283791633342477, e2.coulomb(cgtoS0, cgtoS0, cgtoS0, cgtoS0), delta);
	}

	@Test
	void ss01() {
		assertEquals(0.8427007900292194, e2.coulomb(cgtoS0, cgtoS0, cgtoS1, cgtoS1), delta);
	}

	@Test
	void sp00() {
		assertEquals(0.9403159699467084, e2.coulomb(cgtoS0, cgtoS0, cgtoP0, cgtoP0), delta);
	}

	@Test
	void sd00() {
		assertEquals(0.8086717345109515, e2.coulomb(cgtoS0, cgtoS0, cgtoD0, cgtoD0), delta);
	}

	@Test
	void sf00() {
		assertEquals(0.7132968291998493, e2.coulomb(cgtoS0, cgtoS0, cgtoF0, cgtoF0), delta);
	}

	@Test
	void pd00() {
		assertEquals(0.8583741496984647, e2.coulomb(cgtoP0, cgtoP0, cgtoD0, cgtoD0), delta);
	}

	@Test
	void dd00() {
		assertEquals(0.854418854766963, e2.coulomb(cgtoD0, cgtoD0, cgtoD0, cgtoD0), delta);
	}

	@Test
	void ff00() {
		assertEquals(0.8186960564969021, e2.coulomb(cgtoF0, cgtoF0, cgtoF0, cgtoF0), delta);
	}

	@Test
	void vrrS0000() {
		assertEquals(4.37335456733,
				e2.vrr(new Vector3D(0, 0, 0), 1.0, new Power(0, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(0, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0, 0),
				delta);
	}

	@Test
	void vrrP0000() {
		assertEquals(0.182223107579,
				e2.vrr(new Vector3D(0, 0, 0), 1.0, new Power(1, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(1, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0, 0),
				delta);
	}

	@Test
	void vrrD0000() {
		assertEquals(0.223223306785,
				e2.vrr(new Vector3D(0, 0, 0), 1.0, new Power(2, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(2, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0, 0),
				delta);
	}

	@Test
	void vrrP1000() {
		assertEquals(-5.63387712455e-06,
				e2.vrr(new Vector3D(1, 2, 3), 1.0, new Power(1, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(1, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0, 0),
				delta);
	}

	@Test
	void vrrD1000() {
		assertEquals(0.00022503308545040895,
				e2.vrr(new Vector3D(1, 2, 3), 1.0, new Power(2, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(2, 0, 0), 1.0, new Vector3D(0, 0, 0), 1.0, 1.0, 0),
				delta);
	}

	// --- coulombRepulsion with non-zero bPower / dPower ---
	//
	// coulomb(ContractedGaussian) always passes through contractedHrr, which
	// reduces b and d to s-type before calling vrr — so it was never broken.
	// coulombRepulsion(primitive params) is called directly during analytic
	// gradient evaluation.

	/**
	 * (s, p_x | s, s) with all centres at the origin equals zero by symmetry.
	 * Before the fix, coulombRepulsion dropped bPower and returned the (s,s|s,s)
	 * value (~4.37) instead.
	 */
	@Test
	void coulombRepulsionBPowerXIsZeroAtOrigin() {
		assertEquals(0.0,
				e2.coulombRepulsion(
						new Vector3D(0, 0, 0), 1.0, new Power(0, 0, 0), 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(1, 0, 0), 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(0, 0, 0), 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(0, 0, 0), 1.0),
				delta);
	}

	/**
	 * (s, s | s, p_x) with all centres at the origin equals zero by symmetry.
	 */
	@Test
	void coulombRepulsionDPowerXIsZeroAtOrigin() {
		assertEquals(0.0,
				e2.coulombRepulsion(
						new Vector3D(0, 0, 0), 1.0, new Power(0, 0, 0), 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(0, 0, 0), 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(0, 0, 0), 1.0,
						new Vector3D(0, 0, 0), 1.0, new Power(1, 0, 0), 1.0),
				delta);
	}

	/**
	 * Verifies the HRR algebraic identity for bPower:
	 *   (a, b_x | c, d) = (a+x, b_x-1 | c, d) + (Ax - Bx)(a, b_x-1 | c, d)
	 * Self-contained — uses no external reference oracle.
	 * With a=(1,2,3) and b=origin, Ax - Bx = 1.
	 */
	@Test
	void coulombRepulsionSatisfiesHrrIdentityForBPower() {
		Vector3D a = new Vector3D(1, 2, 3);
		Vector3D b = new Vector3D(0, 0, 0);
		Vector3D c = new Vector3D(0, 0, 0);
		Vector3D d = new Vector3D(0, 0, 0);
		double n = 1.0, alpha = 1.0;

		// LHS: b has Power(1,0,0) — the broken case
		double lhs = e2.coulombRepulsion(
				a, n, new Power(0, 0, 0), alpha,
				b, n, new Power(1, 0, 0), alpha,
				c, n, new Power(0, 0, 0), alpha,
				d, n, new Power(0, 0, 0), alpha);

		// RHS: HRR expansion — bPower reduced to zero, no HRR needed
		double rhs = e2.coulombRepulsion(a, n, new Power(1, 0, 0), alpha,
				b, n, new Power(0, 0, 0), alpha, c, n, new Power(0, 0, 0), alpha,
				d, n, new Power(0, 0, 0), alpha)
				+ (a.getX() - b.getX()) * e2.coulombRepulsion(a, n, new Power(0, 0, 0), alpha,
						b, n, new Power(0, 0, 0), alpha, c, n, new Power(0, 0, 0), alpha,
						d, n, new Power(0, 0, 0), alpha);

		assertEquals(lhs, rhs, delta);
	}

}
