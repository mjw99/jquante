package name.mjw.jquante.common.resource;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
/**
 * 
 * @author mw529
 * Example here: http://www.vogella.com/tutorials/JAXB/article.html
 */
public class NWChemBasisSetFileWriter {

	public static void main(String[] args) throws JAXBException {

		NWChemBasisSetFile nWChemBasisSetFile = new NWChemBasisSetFile();

		nWChemBasisSetFile
				.read("/usr/local/shared/ubuntu-16.04/x86_64/nwchem/6.6/data/libraries/aug-cc-pvtz");
		
		//nWChemBasisSetFile.read("/tmp/aug-cc-pvtz");

		// create JAXB context and instantiate marshaller
		JAXBContext context = null;
		context = JAXBContext.newInstance(NWChemBasisSetFile.class);

		Marshaller m = context.createMarshaller();
		m.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

		// Write to System.out
		m.marshal(nWChemBasisSetFile, System.out);

	}

}
