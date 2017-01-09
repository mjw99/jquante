package name.mjw.jquante.math.optimizer.impl;

import name.mjw.jquante.math.optimizer.FunctionPoint1D;

/**
 * A set of three points a, b, c that bracket a minimum of a function f(x). The
 * following relationships are guaranteed:
 * 
 * <ul>
 * <li>a &lt; b &lt; c</li>
 * <li>f(b) &lt; f(a)</li>
 * <li>f(b) &lt; f(c)</li>
 * </ul>
 * 
 * 
 * @author J. Milthorpe
 * @version 2.0 (Part of MeTA v2.0)
 */
public class Bracket {
	public FunctionPoint1D a, b, c;

	public Bracket(FunctionPoint1D a, FunctionPoint1D b, FunctionPoint1D c) {
		this.a = a;
		this.b = b;
		this.c = c;
	}
}
