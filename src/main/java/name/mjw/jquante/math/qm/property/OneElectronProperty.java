package name.mjw.jquante.math.qm.property;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import name.mjw.jquante.common.Utility;
import name.mjw.jquante.math.qm.SCFMethod;
import name.mjw.jquante.molecule.property.electronic.GridProperty;
import name.mjw.jquante.molecule.property.electronic.PointProperty;
import name.mjw.jquante.parallel.AbstractSimpleParallelTask;
import name.mjw.jquante.parallel.SimpleParallelTask;
import name.mjw.jquante.parallel.SimpleParallelTaskExecuter;

/**
 * Abstract representation of an one-electron property.
 * 
 * @author V.Ganesh
 * @version 2.0 (Part of MeTA v2.0)
 */
public abstract class OneElectronProperty {

	/**
	 * Creates a new instance of OneElectronProperty
	 * */
	public OneElectronProperty() {
		this.scfMethod = null;
	}

	/**
	 * Creates a new instance of OneElectronProperty.
	 *
	 * @param scfMethod
	 *            the Self Consistent Field (SCF) method.
	 */
	public OneElectronProperty(SCFMethod scfMethod) {
		this.scfMethod = scfMethod;
	}

	protected SCFMethod scfMethod;

	/**
	 * Get the value of scfMethod
	 * 
	 * @return the value of scfMethod
	 */
	public SCFMethod getScfMethod() {
		return scfMethod;
	}

	/**
	 * Set the value of scfMethod
	 *
	 * @param scfMethod
	 *            new value of scfMethod
	 */
	public void setScfMethod(SCFMethod scfMethod) {
		this.scfMethod = scfMethod;
	}

	/**
	 * Computes the one electron property on the specified point and returns its
	 * value at the specified point. <br>
	 * Note that the unit of Point3D object must be a.u. No attempt is made to
	 * verify this.
	 * 
	 * @param point
	 *            the point of interest
	 * @return the value of this property at this point
	 */
	public abstract double compute(Vector3D point);

	/**
	 * Computes the one electron property on the specified points and returns
	 * its value at the specified points. <br>
	 * Note that the unit of Point3D object must be a.u. No attempt is made to
	 * verify this.
	 *
	 * @param points
	 *            the points of interest
	 * @return the values of this property at each of the point
	 */
	public double[] compute(Vector3D[] points) {
		double[] fValues = new double[points.length];

		SimpleParallelTask spt = new ComputeOverPoints(points, fValues);
		SimpleParallelTaskExecuter spte = new SimpleParallelTaskExecuter();

		spte.execute(spt);

		return fValues;
	}

	/** Parallel task description to compute the property over a set of points */
	private class ComputeOverPoints extends AbstractSimpleParallelTask {

		private Vector3D[] points;
		private double[] fValues;

		public ComputeOverPoints(Vector3D[] points, double[] fValues) {
			this.points = points;
			this.fValues = fValues;

			setTaskName("OneElectron property computation thread for: "
					+ OneElectronProperty.class.toString());
			setTotalItems(points.length);
		}

		@Override
		public SimpleParallelTask init(int startItem, int endItem) {
			ComputeOverPoints cop = new ComputeOverPoints(points, fValues);

			cop.startItem = startItem;
			cop.endItem = endItem;

			return cop;
		}

		@Override
		public void run() {
			for (int i = startItem; i < endItem; i++)
				fValues[i] = compute(points[i]);
		}
	}

	/**
	 * Computes the one-electron property specified by the GridProperty object. <br>
	 * The units of GridProperty are automatically converted to a.u. It is
	 * always assumed that GridProperty is specified in angstroms.
	 * 
	 * @param gp
	 *            the GridProperty object describing the region of interest to
	 *            compute the one-electron properties
	 */
	public void compute(GridProperty gp) {
		int nx = gp.getNoOfPointsAlongX(), ny = gp.getNoOfPointsAlongY(), nz = gp
				.getNoOfPointsAlongZ();
		double xinc = gp.getXIncrement(), yinc = gp.getYIncrement(), zinc = gp
				.getZIncrement();

		Vector3D ul = gp.getBoundingBox().getUpperLeft();
		double xmin = ul.getX(), ymin = ul.getY(), zmin = ul.getZ();

		Vector3D[] points = new Vector3D[nx * ny * nz];

		double x, y, z;
		int i, j, k, ii;

		ii = 0;
		for (i = 0; i < nx; i++) {
			x = (xmin + (i * xinc)) / Utility.AU_TO_ANGSTROM_FACTOR;
			for (j = 0; j < ny; j++) {
				y = (ymin + (j * yinc)) / Utility.AU_TO_ANGSTROM_FACTOR;
				for (k = 0; k < nz; k++) {
					z = (zmin + (k * zinc)) / Utility.AU_TO_ANGSTROM_FACTOR;

					points[ii++] = new Vector3D(x, y, z);
				}
			}
		}

		gp.setFunctionValues(compute(points));
	}

	/**
	 * Computes the one-electron property specified by the PointProperty object. <br>
	 * Note that the unit of Point3D object must be a.u. No attempt is made to
	 * verify this.
	 * 
	 * @param pp
	 *            the PointProperty object describing the region of interest to
	 *            compute the one-electron properties
	 */
	public void compute(PointProperty pp) {
		pp.setValue(compute(pp.getPoint()));
	}
}
