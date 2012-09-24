package name.mjw.jquante;

public class Atom {

	// Atom id?
	private int atid;

	// Z number
	private int atno;

	// Position
	private float[] r = new float[3];

	// Force
	private float[] f = new float[3];

	// Velocity
	private float[] v = new float[3];

	public Atom(int atno, float x, float y, float z, int atid, float fx,
			float fy, float fz, float vx, float vy, float vz) {
		this.atno = atno;

		this.r[0] = x;
		this.r[1] = y;
		this.r[2] = z;

		this.atid = atid;

		this.f[0] = fx;
		this.f[1] = fy;
		this.f[2] = fz;

		this.v[0] = vx;
		this.v[1] = vy;
		this.v[2] = vz;

	}

	// Set Cartesian position
	public float[] getPos() {
		return r;
	}

	// Get Cartesian position
	public void setPos(float[] r) {
		this.r = r;
	}

	// Get force
	public float[] getForce() {
		return f;
	}

	// Set force
	public void setForce(float[] f) {
		this.f = f;
	}

	// Get velocity
	public float[] getVelocity() {
		return v;
	}

	// Set velocity
	public void setVelocity(float[] v) {
		this.v = v;
	}

	@Override
	public String toString() {
		StringBuilder result = new StringBuilder();
		String NEW_LINE = System.getProperty("line.separator");

		result.append(String.format(
				" Atom ID: %d Atomic Num: %2d (%6.3f,%6.3f,%6.3f)", atid, atno,
				r[0], r[1], r[2]));

		result.append(NEW_LINE);

		return result.toString();

	}

}
