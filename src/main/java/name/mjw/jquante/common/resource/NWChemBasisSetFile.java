package name.mjw.jquante.common.resource;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Locale;

import java.util.List;

import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.log4j.Logger;

import name.mjw.fortranformat.FortranFormat;

@XmlRootElement(name = "basis")
public class NWChemBasisSetFile {
	private final static Logger LOG = Logger
			.getLogger(NWChemBasisSetFile.class);

	@XmlAttribute(name = "name")
	private String basisSetName;

	@XmlElement(name = "atom")
	private List<Library> atoms;

	/**
	 * Representation of a NWChem basis file.
	 * 
	 * 
	 * @see <a
	 *      href="https://github.com/jeffhammond/nwchem/blob/HEAD/src/basis/basis_dox.c">src/basis/basis_dox.c</a>
	 * 
	 * 
	 */
	public NWChemBasisSetFile() {
		atoms = new ArrayList<Library>();
	}

	/**
	 * Reads a NWChem basis file.
	 * 
	 * @param fileName
	 *            Filename of NWChem basis file.
	 * @throws ParseException
	 * @throws IOException
	 * 
	 * @see <a
	 *      href="https://github.com/jeffhammond/nwchem/blob/HEAD/src/basis/bas_input.F">src/basis/bas_input.F</a>
	 * 
	 * 
	 */
	public void read(String fileName) {

		basisSetName = fileName;

		try (BufferedReader in = new BufferedReader(new FileReader(fileName))) {

			String line;

			while ((line = in.readLine()) != null) {

				if (line.matches("^basis.*")) {

					parseBasisSection(line, in);

				}

			}

			in.close();

		} catch (IOException e) {
			e.printStackTrace();

		} catch (ParseException e) {
			e.printStackTrace();
		}
	}

	/**
	 * Parse the section:
	 * 
	 * <pre>
	 * basis "He_aug-cc-pVTZ" SPHERICAL 
	 * ... 
	 * end
	 * </pre>
	 * 
	 * 
	 * @param line
	 * @param br
	 * @throws ParseException
	 * @throws IOException
	 */
	void parseBasisSection(String line, BufferedReader br)
			throws ParseException, IOException {

		String elementName;
		String basisName;

		String words[];

		words = line.split("\\s+");

		// basis "Na_cc-pVTZ" SPHERICAL
		LOG.debug(words[0] + " " + words[1] + " " + words[2]);

		if (!(words[0].contentEquals("basis"))) {
			throw new ParseException("Expecting basis term", 1);
		}

		elementName = words[1].replace("\"", "").split("_")[0];
		basisName = words[1];

		List<Shell> orbitals = new ArrayList<Shell>();

		// while ((br.readLine()) != null) {
		//
		// if (!(line.matches("^end.*"))) {
		//
		// orbitals.addAll(parseOrbitalSection(br));
		//
		// }
		//
		// }

		orbitals.addAll(parseShellSection(br));

		Library atom = new Library(elementName);
		atom.addOrbitals(orbitals);
		atoms.add(atom);
	}

	/**
	 * Parse the section:
	 * 
	 * <pre>
	 * H    S
	 *      33.8700000              0.0060680
	 *       5.0950000              0.0453080
	 *       1.1590000              0.2028220
	 * </pre>
	 * 
	 * 
	 * @param line
	 * @param br
	 * @throws ParseException
	 * @throws IOException
	 */
	static List<Shell> parseShellSection(BufferedReader br)
			throws ParseException, IOException {

		int numberOfShellsToParse = 1;
		Double exp;
		Double coeff;
		String line;
		ArrayList<Shell> shells = null;
		ArrayList<Object> lineObjects;
		final String orbitalFormat = "(1X, E18.7, 9X, 6(E15.7, 8X))";

		DecimalFormat df = new DecimalFormat("0",
				DecimalFormatSymbols.getInstance(Locale.ENGLISH));
		df.setMaximumFractionDigits(340); // 340 =
											// DecimalFormat.DOUBLE_FRACTION_DIGITS

		line = br.readLine();
		lineObjects = FortranFormat.read(line, "(2A5)");

		String currentElement = (String) lineObjects.get(0);
		String shellType = (String) lineObjects.get(1);
		//TODO check for a SP shell and ensure 

		// Peek ahead here to work out how many Shells to return
		br.mark(5000);
		line = br.readLine();

		// Hack since FortranFormat does not
		// understand Double-Precision Exponent
		line = line.replace("D", "E");
		lineObjects = FortranFormat.read(line, orbitalFormat);
		numberOfShellsToParse = getNonNullLength(lineObjects) - 1;

		br.reset();

		shells = new ArrayList<Shell>(numberOfShellsToParse);

		for (int i = 0; i < numberOfShellsToParse; i++) {
			shells.add(new Shell(shellType));
		}

		while (true) {

			br.mark(5000);
			line = br.readLine();

			/**
			 * Time to switch to a new Shell(shellType)
			 */
			if (line.contains(currentElement)) {
				br.reset();
				break;
			}

			/**
			 * Time to switch to a new
			 */
			if (line.contains("end")) {
				break;
			}

			// Hack since FortranFormat does not understand Double-Precision
			// Exponent
			line = line.replace("D", "E");

			lineObjects = FortranFormat.read(line, orbitalFormat);

			exp = (Double) lineObjects.get(0);

			for (int i = 0; i < numberOfShellsToParse; i++) {

				coeff = (Double) lineObjects.get(i + 1);

				shells.get(i).addExponentCoefficientPair(df.format(exp),
						df.format(coeff));
			}
		}

		return shells;
	}

	static class Library {
		@XmlAttribute(name = "symbol")
		String element;

		@XmlAttribute
		int atomicNumber;

		@XmlElement(name = "orbital")
		List<Shell> orbitals;

		Library(String element) {
			this.element = element;
			orbitals = new ArrayList<Shell>();
		}

		void addOrbitals(List<Shell> orbitals) {
			this.orbitals.addAll(orbitals);
		}

	}

	static class Shell {

		@XmlAttribute(name = "type")
		String type;

		@XmlElement(name = "entry")
		List<exponentCoefficientPair> exponentCoefficientPairs = null;

		Shell(String type) {
			this.type = type;
			exponentCoefficientPairs = new ArrayList<exponentCoefficientPair>();
		}

		void addExponentCoefficientPair(String exponent, String coefficient) {

			exponentCoefficientPair pair = new exponentCoefficientPair(
					exponent, coefficient);

			exponentCoefficientPairs.add(pair);

		}

	}

	static class exponentCoefficientPair {

		@XmlAttribute(name = "exp")
		String exponent;

		@XmlAttribute(name = "coeff")
		String coefficient;

		exponentCoefficientPair(String exponent, String coefficient) {
			this.exponent = exponent;
			this.coefficient = coefficient;
		}

	}

	/**
	 * Count the number of non-null elements in an array
	 * 
	 * @param lineObjects
	 * @return number of non-null elements
	 */
	private static int getNonNullLength(ArrayList<Object> lineObjects) {
		int count = 0;
		for (Object el : lineObjects)
			if (el != null)
				++count;
		return count;
	}
}
