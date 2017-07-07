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
			.getLogger(NWChemToXMLConverter.class);

	@XmlAttribute(name = "name")
	private String basisSetName;

	@XmlElement(name = "atom")
	private List<Atom> atoms;

	public NWChemBasisSetFile() {
		atoms = new ArrayList<Atom>();
	}

	/**
	 * Reads a NWChem basis file.
	 * 
	 * @param fileName
	 *            Filename of NWChem basis file.
	 * @throws ParseException
	 * @throws IOException
	 */
	public void read(String fileName) {

		basisSetName = fileName;

		try {
			BufferedReader in = new BufferedReader(new FileReader(fileName));

			String line;

			while ((line = in.readLine()) != null) {

				if (line.matches("^basis.*")) {

					parseAtomSection(line, in);

				}

			}

			in.close();

		} catch (IOException e) {
			e.printStackTrace();

		} catch (ParseException e) {
			e.printStackTrace();
		}
	}

	void parseAtomSection(String line, BufferedReader br)
			throws ParseException, IOException {

		String elementName;
		String elementBasisSetName;

		String words[];

		words = line.split("\\s+");

		// basis "Na_cc-pVTZ" SPHERICAL
		LOG.debug(words[0] + " " + words[1] + " " + words[2]);

		if (!(words[0].contentEquals("basis"))) {
			throw new ParseException("Expecting basis term", 1);
		}

		elementName = words[1].replace("\"", "").split("_")[0];
		elementBasisSetName = words[1];

		List<Orbital> orbitals = new ArrayList<Orbital>();
		orbitals.addAll(parseOrbitalSection(br));

		Atom atom = new Atom(elementName);
		atom.addOrbitals(orbitals);
		atoms.add(atom);

	}

	static List<Orbital> parseOrbitalSection(BufferedReader br)
			throws ParseException, IOException {

		int numberOfOrbitalsToParse = 1;
		Double exp;
		Double coeff;
		String line;
		ArrayList<Orbital> orbitals = null;
		ArrayList<Object> lineObjects;
		final String orbitalFormat = "(1X, E18.7, 9X, 6(E15.7, 8X))";
		

		DecimalFormat df = new DecimalFormat("0",
				DecimalFormatSymbols.getInstance(Locale.ENGLISH));
		df.setMaximumFractionDigits(340); // 340 =
											// DecimalFormat.DOUBLE_FRACTION_DIGITS

		line = br.readLine();
		lineObjects = FortranFormat.read(line, "(2A5)");

		String currentElement = (String) lineObjects.get(0);
		String orbitalType = (String) lineObjects.get(1);

		// Peek ahead here to work out how many orbitalElement to return now
		br.mark(5000);
		line = br.readLine();

		// Hack since FortranFormat does not understand Double-Precision
		// Exponent
		line = line.replace("D", "E");
		lineObjects = FortranFormat.read(line, orbitalFormat);
		numberOfOrbitalsToParse = getNonNullLength(lineObjects) - 1;

		br.reset();

		orbitals = new ArrayList<Orbital>(numberOfOrbitalsToParse);

		for (int i = 0; i < numberOfOrbitalsToParse; i++) {
			orbitals.add(new Orbital(orbitalType));
		}

		while (true) {

			line = br.readLine();

			/**
			 * Time to switch to a new Orbital(orbitalType)
			 */
			if (line.contains(currentElement)) {
				break;

//			System.out.println(line);
//				br.reset();
//				parseOrbitalSection(br);
				//return orbitals;
//				lineObjects = FortranFormat.read(line, "(2A5)");
//				orbitalType = (String) lineObjects.get(1);
//				
//				br.mark(5000);
//				line = br.readLine();
//				
//				
//				line = line.replace("D", "E");
//				System.out.println(line);
//				lineObjects = FortranFormat.read(line, orbitalFormat);
//				numberOfOrbitalsToParse = getNonNullLength(lineObjects) - 1;
//				
//				//Add extra orbitals
//				for (int i = 0; i < numberOfOrbitalsToParse; i++) {
//					orbitals.add(new Orbital(orbitalType));
//				}
//								continue;
				

				
			}

			if (line.contains("end")) {
				break;				
			}

			// Hack since FortranFormat does not understand Double-Precision
			// Exponent
			line = line.replace("D", "E");

			lineObjects = FortranFormat.read(line, orbitalFormat);

			exp = (Double) lineObjects.get(0);

			for (int i = 0; i < numberOfOrbitalsToParse; i++) {

				coeff = (Double) lineObjects.get(i + 1);

				orbitals.get(i).addExponentCoefficientPair(df.format(exp),
						df.format(coeff));
			}
		}

		return orbitals;
	}

	static class Atom {
		@XmlAttribute(name = "symbol")
		String element;

		@XmlAttribute
		int atomicNumber;

		@XmlElement(name = "orbital")
		List<Orbital> orbitals;

		Atom(String element) {
			this.element = element;
			orbitals = new ArrayList<Orbital>();
		}

		void addOrbitals(List<Orbital> orbitals) {

			this.orbitals.addAll(orbitals);

		}

	}

	static class Orbital {

		@XmlAttribute(name = "type")
		String type;

		@XmlElement(name = "member")
		List<exponentCoefficientPair> exponentCoefficientPairs = null;

		Orbital(String type) {
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
	public static int getNonNullLength(ArrayList<Object> lineObjects) {
		int count = 0;
		for (Object el : lineObjects)
			if (el != null)
				++count;
		return count;
	}
}
