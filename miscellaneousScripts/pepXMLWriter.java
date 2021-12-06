package Features;

import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsPipelineAnalysis;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import java.nio.file.Path;
import java.nio.file.Paths;

public class pepXMLWriter {

    //probably should be put in new class pepxmlWriter
    public static void writepepXML(MsmsPipelineAnalysis newMsmsPipelineAnalysis, String outFile) throws JAXBException {
        //write xml
        JAXBContext jaxbContext = JAXBContext.newInstance(MsmsPipelineAnalysis.class);
        Marshaller marshaller = jaxbContext.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

        Path path = Paths.get(outFile);
        java.io.File file = path.toFile();
        marshaller.marshal(newMsmsPipelineAnalysis, file);
    }
}
