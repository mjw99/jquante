package name.mjw.jquante.common.resource;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.nio.file.Files;
import java.nio.file.Paths;

import java.security.DigestInputStream;
import java.security.MessageDigest;
import java.security.NoSuchAlgorithmException;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Locale;

import java.util.List;

import javax.xml.bind.DatatypeConverter;
import javax.xml.bind.annotation.XmlAttribute;
import javax.xml.bind.annotation.XmlElement;
import javax.xml.bind.annotation.XmlRootElement;

import org.apache.log4j.Logger;

import name.mjw.fortranformat.FortranFormat;
import name.mjw.jquante.config.impl.AtomInfo;

@XmlRootElement(name = "basis")
public class NWChemBasisSetFile {
	private final static Logger LOG = Logger
			.getLogger(NWChemBasisSetFile.class);

	private final static int markLength = 5000;

	@XmlAttribute(name = "name")
	private String basisSetName;

	@XmlAttribute(name = "md5sum")
	private String md5sumOfInputBasisFile;

	@XmlElement(name = "atom")
	private List<Library> atoms;

	/**
	 * A representation of a NWChem basis file.
	 * 
	 * 
	 * @see <a
	 *      href="https://github.com/jeffhammond/nwchem/blob/HEAD/src/basis/basis_dox.c">src/basis/basis_dox.c</a>
	 * @see <a
	 *      href="http://www.ccl.net/cca/documents/basis-sets/basis.html">General
	 *      basis set naming information</a>
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
	 * @throws NoSuchAlgorithmException
	 * @throws ParseException
	 * @throws IOException
	 * 
	 * @see <a
	 *      href="https://github.com/jeffhammond/nwchem/blob/HEAD/src/basis/bas_input.F">src/basis/bas_input.F</a>
	 * 
	 * 
	 */
	public void read(String fileName) {

		File f = new File(fileName);
		basisSetName = f.getName();
		MessageDigest md = null;

		try {
			md = MessageDigest.getInstance("MD5");
		} catch (NoSuchAlgorithmException e1) {
			e1.printStackTrace();
		}

		// https://stackoverflow.com/questions/304268/getting-a-files-md5-checksum-in-java
		try (InputStream is = Files.newInputStream(Paths.get(fileName));
				DigestInputStream dis = new DigestInputStream(is, md);
				BufferedReader in = new BufferedReader(new InputStreamReader(
						dis, "UTF-8"))) {

			String line;

			while ((line = in.readLine()) != null) {

				if (line.matches("^basis.*")) {

					parseBasisSection(line, in);

				}

			}

			byte[] digest = md.digest();
			md5sumOfInputBasisFile = DatatypeConverter.printHexBinary(digest)
					.toLowerCase();

		} catch (IOException | ParseException e) {
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

		Library atom = new Library(elementName);

		List<Shell> shells = new ArrayList<Shell>();

		br.mark(markLength);
		while ((line = br.readLine()) != null) {

			if (line.matches("^end.*")) {

				br.reset();
				break;

			} else {

				LOG.debug("dispatching to parseShellSection() with line ");
				LOG.debug(line);

				br.reset();
				shells.addAll(parseShellSection(br, elementName));

			}

			br.mark(markLength);

		}

		atom.addShell(shells);
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
	static List<Shell> parseShellSection(BufferedReader br, String elementName)
			throws ParseException, IOException {

		int numberOfContractedShellsToParse = 1;
		Double exp;
		Double coeff;
		String line;
		ArrayList<Shell> shells = null;
		ArrayList<Object> lineObjects;
		final String shellFormat = "(1X, E18.7, 9X, 6(E15.7, 8X))";

		DecimalFormat df = new DecimalFormat("0",
				DecimalFormatSymbols.getInstance(Locale.ENGLISH));
		df.setMaximumFractionDigits(340); // 340 =
		// DecimalFormat.DOUBLE_FRACTION_DIGITS

		line = br.readLine();
		LOG.debug("");
		LOG.debug("parseShellSection() is starting on line: ");
		LOG.debug(line);

		lineObjects = FortranFormat.read(line, "(2A5)");

		String currentElement = (String) lineObjects.get(0);
		String shellType = (String) lineObjects.get(1);

		// Peek ahead here to work out how many contractedShells to return
		br.mark(markLength);
		line = br.readLine();

		// Hack since FortranFormat does not
		// understand Double-Precision Exponent
		line = line.replace("D", "E");
		lineObjects = FortranFormat.read(line, shellFormat);
		numberOfContractedShellsToParse = getNonNullLength(lineObjects) - 1;
		br.reset();

		shells = new ArrayList<Shell>(numberOfContractedShellsToParse);

		// For the basis sets in which s- and p-type functions share the same
		// exponents, the term SP-shell is used.
		if (shellType.equals("SP")) {
			LOG.debug("SP Found");

			for (int i = 0; i < numberOfContractedShellsToParse; i++) {
				LOG.debug(ShellType.values()[i]);
				shells.add(new Shell(ShellType.values()[i].toString()));

			}
		} else {
			for (int i = 0; i < numberOfContractedShellsToParse; i++) {
				shells.add(new Shell(shellType));
			}
		}

		while (true) {
			line = br.readLine();

			if (line.matches("^" + elementName + ".*") | line.matches("^end.*")) {

				LOG.debug("Seen line: " + line);
				LOG.debug("New shell time..");
				LOG.debug("");

				br.reset();
				return shells;
			}

			// Hack since FortranFormat does not understand Double-Precision
			// Exponent
			line = line.replace("D", "E");

			lineObjects = FortranFormat.read(line, shellFormat);
			LOG.debug(lineObjects);

			exp = (Double) lineObjects.get(0);

			for (int i = 0; i < numberOfContractedShellsToParse; i++) {

				coeff = (Double) lineObjects.get(i + 1);

				shells.get(i).addExponentCoefficientPair(df.format(exp),
						df.format(coeff));
			}

			br.mark(markLength);
		}

	}

	static class Library {
		@XmlAttribute(name = "symbol")
		String element;

		@XmlAttribute
		int atomicNumber;

		@XmlElement(name = "orbital")
		List<Shell> shells;

		Library(String element) {
			this.element = element;
			atomicNumber = AtomInfo.getInstance().getAtomicNumber(element);
			shells = new ArrayList<Shell>();
		}

		void addShell(List<Shell> orbitals) {
			this.shells.addAll(orbitals);
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

	private enum ShellType {
		S(0), P(1), D(2), F(3), G(4), H(5);

		private int value;

		private ShellType(int value) {
			this.value = value;
		}
	}
}
