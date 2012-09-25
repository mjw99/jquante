package name.mjw.jquante;

public class GaussianTypeOrbital {

	
	/*
	 *  g(x,y,z) = A*(x^i)*(y^j)*(z^k)*exp{-a*(r-ro)^2}
	 */
	
	
	private double alpha;
	
	// Origin
	private double x = 0.0;
	private double y = 0.0;
	private double z = 0.0;
	
	private int l = 0;
	private int m = 0;
	private int n = 0;
	
	private double norm;
	private double coef = 1.0;
	
	public GaussianTypeOrbital(double alpha, double x, double y, double z, int l, int m,
			int n, double coef) {
	
		this.alpha = alpha;

		this.x = x;
		this.y = y;
		this.z = z;

		this.l = l;
		this.m = m;
		this.n = n;

		this.coef = coef;

	}
	
	public void normalise() {
		
	}
	
	public void overlap(GaussianTypeOrbital Other) {
		
	}
	
	public void amp(double x, double y, double z) {
		
	}

	
	
	
}
