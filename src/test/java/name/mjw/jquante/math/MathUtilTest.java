package name.mjw.jquante.math;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealVector;
import org.junit.jupiter.api.Test;

class MathUtilTest {
	double diff = 0.00001;

	@Test
	void factorial2MinusOne() {

		assertEquals(1, MathUtil.factorial2(-1));

	}

	@Test
	void factorial2Zero() {

		assertEquals(1, MathUtil.factorial2(0));

	}

	@Test
	void factorial2Nine() {

		assertEquals(945, MathUtil.factorial2(9));

	}

	@Test
	void testfindDihedralAngle1() {

		Vector3D i = new Vector3D(60.614941, 67.371346, 12.713049);
		Vector3D j = new Vector3D(59.104156, 67.044731, 12.827555);
		Vector3D k = new Vector3D(58.827194, 65.911911, 13.835642);
		Vector3D l = new Vector3D(59.367950, 66.384277, 15.233553);

		assertEquals(59.66383, Math.toDegrees(MathUtil.findDihedral(i, j, k, l)), diff);

	}

	@Test
	void testfindDihedralAngle2() {

		Vector3D i = new Vector3D(30.958572, 73.709137, 38.031029);
		Vector3D j = new Vector3D(31.924915, 72.577698, 37.985279);
		Vector3D k = new Vector3D(31.643818, 71.563255, 36.873047);
		Vector3D l = new Vector3D(31.253489, 70.131218, 37.371773);

		assertEquals(-115.172063, Math.toDegrees(MathUtil.findDihedral(i, j, k, l)), diff);

	}

	@Test
	void findAngleRightAngle() {
		// v1=(1,0,0), v2=(0,0,0), v3=(0,1,0) => angle at v2 = 90 degrees
		Vector3D v1 = new Vector3D(1, 0, 0);
		Vector3D v2 = new Vector3D(0, 0, 0);
		Vector3D v3 = new Vector3D(0, 1, 0);
		assertEquals(90.0, Math.toDegrees(MathUtil.findAngle(v1, v2, v3)), diff);
	}

	@Test
	void findAngleStraight() {
		// v1=(1,0,0), v2=(0,0,0), v3=(-1,0,0) => 180 degrees
		Vector3D v1 = new Vector3D(1, 0, 0);
		Vector3D v2 = new Vector3D(0, 0, 0);
		Vector3D v3 = new Vector3D(-1, 0, 0);
		assertEquals(180.0, Math.toDegrees(MathUtil.findAngle(v1, v2, v3)), diff);
	}

	@Test
	void factorialRatioSquaredKnownValue() {
		// a=2, b=1: 2! / 1! / (2-2*1)! = 2/1/1 = 2.0
		assertEquals(2.0, MathUtil.factorialRatioSquared(2, 1), diff);
	}

	@Test
	void binomialPrefactorKnownValue() {
		// s=0, ia=0, ib=0, xpa=1.0, xpb=1.0
		// loop t=0: condition (0-0)<=0 && 0<=0 => C(0,0)*C(0,0)*1^0*1^0 = 1
		assertEquals(1.0, MathUtil.binomialPrefactor(0, 0, 0, 1.0, 1.0), diff);
	}

	@Test
	void realMatrixToRealVectorRowMajorOrder() {
		Array2DRowRealMatrix m = new Array2DRowRealMatrix(new double[][]{{1, 2}, {3, 4}});
		RealVector v = MathUtil.realMatrixToRealVector(m);
		assertEquals(4, v.getDimension());
		assertEquals(1.0, v.getEntry(0), 1e-10);
		assertEquals(2.0, v.getEntry(1), 1e-10);
		assertEquals(3.0, v.getEntry(2), 1e-10);
		assertEquals(4.0, v.getEntry(3), 1e-10);
	}

	@Test
	void matrixToStringReturnsNonNull() {
		Array2DRowRealMatrix m = new Array2DRowRealMatrix(new double[][]{{1.5, -2.3}, {0.0, 4.0}});
		assertNotNull(MathUtil.matrixToString(m));
	}

}
