package name.mjw.jquante.math.qm.basis;

public class CompactContractedGaussian {
	public final int l;
	public final int m;
	public final int n;

	public final double x;
	public final double y;
	public final double z;

	public final double normalization;
	public final double[] exponents;
	public final double[] coefficients;
	public final double[] primnorms;

	CompactContractedGaussian(int l, int m, int n, double x, double y, double z, double normalization, double[] exponents,
			double[] coefficients, double[] primnorms) {

		this.l = l;
		this.m = m;
		this.n = n;

		this.x = x;
		this.y = y;
		this.z = z;
		
		this.normalization = normalization;

		this.exponents = exponents;
		this.coefficients = coefficients;
		this.primnorms = primnorms;

	}
}
