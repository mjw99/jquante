package name.mjw.jquante.math.qm.integral;


import static org.junit.jupiter.api.Assertions.assertEquals;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.junit.jupiter.api.Test;


class IntegralsUtilTest {

	private final double delta = 0.00001;

	@Test
	void gaussianProductCenterOne() {

		assertEquals(new Vector3D(0.0, 0.0, 0.0),
				IntegralsUtil.gaussianProductCenter(1.0, new Vector3D(0, 0, 0), 1.0, new Vector3D(0, 0, 0)));
	}

	@Test
	void gaussianProductCenterTwo() {

		assertEquals(new Vector3D(0.5, 0.5, 0.5),
				IntegralsUtil.gaussianProductCenter(1.0, new Vector3D(0, 0, 0), 1.0, new Vector3D(1, 1, 1)));
	}

	@Test
	void ijkl2intindex0000() {

		assertEquals(0, IntegralsUtil.ijkl2intindex(0, 0, 0, 0));
	}

	@Test
	void ijkl2intindex1000() {

		assertEquals(1, IntegralsUtil.ijkl2intindex(1, 0, 0, 0));
	}

	@Test
	void ijkl2intindex0100() {

		assertEquals(1, IntegralsUtil.ijkl2intindex(0, 1, 0, 0));
	}

	@Test
	void ijkl2intindex0010() {

		assertEquals(1, IntegralsUtil.ijkl2intindex(0, 0, 1, 0));
	}

	@Test
	void ijkl2intindex0001() {

		assertEquals(1, IntegralsUtil.ijkl2intindex(0, 0, 0, 1));
	}

	@Test
	void ijkl2intindex1001() {

		assertEquals(2, IntegralsUtil.ijkl2intindex(1, 0, 0, 1));
	}

	@Test
	void ijkl2intindex1100() {

		assertEquals(3, IntegralsUtil.ijkl2intindex(1, 1, 0, 0));
	}

	@Test
	void ijkl2intindex0011() {

		assertEquals(3, IntegralsUtil.ijkl2intindex(0, 0, 1, 1));
	}

	@Test
	void computeFGamma() {
		assertEquals(1.0, IntegralsUtil.computeFGamma(0, 0), delta);

	}

	@Test
	void gammaIncompleteOne() {
		assertEquals(1.49365, IntegralsUtil.lowerIncompleteGamma(0.5, 1), delta);

	}

	@Test
	void gammaIncompleteTwo() {
		assertEquals(0.6545103, IntegralsUtil.lowerIncompleteGamma(1.5, 2), delta);

	}

	@Test
	void gammaIncompleteThree() {
		assertEquals(0, IntegralsUtil.lowerIncompleteGamma(2.5, 1E-12), delta);

	}

	@Test
	void logGammaAtOne() {
		// Γ(1) = 0! = 1, so ln(Γ(1)) = 0
		assertEquals(0.0, IntegralsUtil.logGamma(1.0), delta);
	}

	@Test
	void logGammaAtTwo() {
		// Γ(2) = 1! = 1, so ln(Γ(2)) = 0
		assertEquals(0.0, IntegralsUtil.logGamma(2.0), delta);
	}

	@Test
	void logGammaAtHalf() {
		// Γ(0.5) = sqrt(π), so ln(Γ(0.5)) = 0.5*ln(π) ≈ 0.57236
		assertEquals(0.57236, IntegralsUtil.logGamma(0.5), delta);
	}

	@Test
	void logGammaAtThree() {
		// Γ(3) = 2! = 2, so ln(Γ(3)) = ln(2) ≈ 0.69315
		assertEquals(0.69315, IntegralsUtil.logGamma(3.0), delta);
	}

	@Test
	void computeFGammaAtM0X1() {
		// F_0(1) = 0.5 * x^{-0.5} * γ(0.5, 1) = 0.5 * 1.0 * 1.49365 ≈ 0.74682
		assertEquals(0.74682, IntegralsUtil.computeFGamma(0, 1.0), delta);
	}

	@Test
	void computeFGammaM0AtZeroApproachesOne() {
		// F_0(0) = 1 (by continuity, as verified by series expansion)
		assertEquals(1.0, IntegralsUtil.computeFGamma(0, 0.0), delta);
	}

	@Test
	void ijkl2intindexSymmetry1111() {
		// (1,1,1,1) with all four equal should map to a unique index
		int idx = IntegralsUtil.ijkl2intindex(1, 1, 1, 1);
		// ij=1*(1+1)/2+1=2, kl=2, ij>=kl, result = 2*3/2+2 = 5
		assertEquals(5, idx);
	}

	@Test
	void ijkl2intindexAllPermutationsAreEqual() {
		// (i,j,k,l) and (k,l,i,j) must produce the same index
		int ij0kl = IntegralsUtil.ijkl2intindex(1, 0, 1, 1);
		int kl0ij = IntegralsUtil.ijkl2intindex(1, 1, 1, 0);
		assertEquals(ij0kl, kl0ij);
	}

}
