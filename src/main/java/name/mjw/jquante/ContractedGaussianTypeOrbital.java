package name.mjw.jquante;

import java.util.HashMap;

import javax.vecmath.Point3d;

public class ContractedGaussianTypeOrbital {

	HashMap<Double, PrimitiveGaussianTypeOrbital> primatives = new HashMap<Double, PrimitiveGaussianTypeOrbital>();

	public void addPrimitiveGaussianTypeOrbital(Double coeff,
			PrimitiveGaussianTypeOrbital prim) {

		primatives.put(coeff, prim);

	}

	/**
	 * Compute the amplitude of the primitive Gaussian at point x,y,z
	 * 
	 * @param point
	 */
	public double valueAtPoint(Point3d point) {
		double value = 0.0;

		for (Double coeff : primatives.keySet()) {
			value += coeff * (primatives.get(coeff).valueAtPoint(point));

		}
		return value;

	}

}
