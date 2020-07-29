import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsRunSummary;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBElement;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Unmarshaller;
import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import javax.xml.transform.stream.StreamSource;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;

public class testRun {

    private static boolean advanceReaderToNextRunSummary(XMLStreamReader xsr)
            throws XMLStreamException {
        do {
            if (xsr.next() == XMLStreamConstants.END_DOCUMENT)
                return false;
        } while (!(xsr.isStartElement() && xsr.getLocalName().equals("msms_run_summary")));

        return true;
    }

    public static void main(String[] args) {
        String file = "23aug2017_hela_serum_timecourse_4mz_narrow_1_rank1.pepXML";
        Path path = Paths.get(file);

        try (FileInputStream fis = new FileInputStream(file)) {
            // we'll manually iterate over msmsRunSummaries - won't need so much memory
            // at once for processing large files.
            JAXBContext ctx = JAXBContext.newInstance(MsmsRunSummary.class);
            Unmarshaller unmarshaller = ctx.createUnmarshaller();

            XMLInputFactory xif = XMLInputFactory.newFactory();

            StreamSource ss = new StreamSource(fis);
            XMLStreamReader xsr = xif.createXMLStreamReader(ss);


            while (advanceReaderToNextRunSummary(xsr)) {
                // we've advanced to the next MsmsRunSummary in the file
                long timeLo = System.nanoTime();
                JAXBElement<MsmsRunSummary> unmarshalled = unmarshaller
                        .unmarshal(xsr, MsmsRunSummary.class);
                long timeHi = System.nanoTime();
                System.out.printf("Unmarshalling took %.4fms (%.2fs)\n",
                        (timeHi-timeLo)/1e6, (timeHi-timeLo)/1e9);
                MsmsRunSummary runSummary = unmarshalled.getValue();
                if (runSummary.getSpectrumQuery().isEmpty()) {
                    String msg = String.format("Parsed msms_run_summary was empty for " +
                                    "'%s', summary base_name '%s'",
                            path.toUri().toString(), runSummary.getBaseName());
                    System.out.println(msg);
                }
            }
        } catch (IOException | JAXBException | XMLStreamException e) {
            e.printStackTrace();
        }
    }
}
