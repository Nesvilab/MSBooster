/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package writers;

import jakarta.xml.bind.JAXBContext;
import jakarta.xml.bind.JAXBException;
import jakarta.xml.bind.Marshaller;
import org.apache.commons.lang3.ArrayUtils;
import peptideptmformatting.PeptideFormatter;
import readers.datareaders.PinReader;
import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.ModificationInfo;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsPipelineAnalysis;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsRunSummary;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SearchHit;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SpectrumQuery;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;

import static utils.Print.printError;
import static utils.Print.printInfo;

public class PepXMLDivider {
    int divisions;
    MsmsRunSummary runSummary;
    MsmsPipelineAnalysis[] analyses;
    ArrayList<SpectrumQuery>[] sqArrays;
    ArrayList<SpectrumQuery>[] sqToBin;
    String pepxmlName;
    HashMap<String, String> pinMap = new HashMap<>();
    HashMap<String, Integer> baseChargetoUnicodeBin = new HashMap<>();

    public PepXMLDivider(int divisions) {
        if (divisions < 2) {
            printError("Divisions for PepXMLDivider must be 2 or greater");
            System.exit(1);
        }
        this.divisions = divisions;
        analyses = new MsmsPipelineAnalysis[divisions];
        sqArrays = new ArrayList[divisions];
        sqToBin = new ArrayList[divisions];
        for (int i = 0; i < divisions; i++) {
            sqArrays[i] = new ArrayList<>();
            sqToBin[i] = new ArrayList<>();
        }
    }

    private void load(String pepxmlPath) throws FileParsingException {
        clear();
        pepxmlName = pepxmlPath;
        Path path = Paths.get(pepxmlPath);
        MsmsPipelineAnalysis analysis = PepXmlParser.parse(path);
        runSummary = analysis.getMsmsRunSummary().get(0);

        //initializing pepxml files
        for (int i = 0; i < divisions; i++) {
            analyses[i] = PepXmlParser.parse(path); //redo so each analysis is a separate entity
            analyses[i].getMsmsRunSummary().get(0).getSpectrumQuery().clear();
        }
    }

    public void process(String pepxmlPath) throws FileParsingException, JAXBException {
        load(pepxmlPath);
        divide();
        writePepXML();
        //can do mapPepxmlToBin
    }

    public void divide() {
        List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();
        for (SpectrumQuery sq : spectrumQueries) {
            int charge = sq.getAssumedCharge();
            SearchHit searchHit = sq.getSearchResult().get(0).getSearchHit().get(0);
            int pepUnicodeSum = charge;
            int bin;

            //get number to identify peptidoform. No unique, but same peptidoform always gets same number
            ModificationInfo mi = searchHit.getModificationInfo();
            String peptide;
            if (mi != null) {
                peptide = mi.getModifiedPeptide();
            } else {
                peptide = searchHit.getPeptide();
            }
            pepUnicodeSum += unicodeSum(peptide);

            bin = pepUnicodeSum % divisions;

            //add sq to every bin except that bin
            for (int i = 0; i < sqArrays.length; i++) {
                if (i != bin) {
                    sqArrays[i].add(sq);
                }
            }

            //opposite of sqArrays
            sqToBin[bin].add(sq);
        }
    }

    private int unicodeSum(String string) {
        int sum = 0;
        char[] chars = string.toCharArray();
        for (char c : chars) {
            sum += c;
        }
        return sum;
    }

    private void writePepXML() throws JAXBException {
        //get file base name
        String[] pathSplit = pepxmlName.split("\\.");
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < pathSplit.length - 1; i++) {
            sb.append(pathSplit[i]).append(".");
        }
        String baseName = sb.toString();

        for (int i = 0; i < divisions; i++) {
            //add spectrum queries
            for (SpectrumQuery sq : sqArrays[i]) {
                analyses[i].getMsmsRunSummary().get(0).getSpectrumQuery().add(sq);
            }

            //write to file (divide{i})
            JAXBContext jaxbContext = JAXBContext.newInstance(MsmsPipelineAnalysis.class);
            Marshaller marshaller = jaxbContext.createMarshaller();
            marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);
            Path path = Paths.get(baseName + "divideMSBooster" + i + ".pepXML");
            java.io.File file = path.toFile();
            marshaller.marshal(analyses[i], file);
        }
    }

    private void clear() {
        for (ArrayList<SpectrumQuery> a : sqArrays) {
            a.clear();
        }
        for (ArrayList<SpectrumQuery> a : sqToBin) {
            a.clear();
        }
    }

    public void loadPin(String pinFile) throws IOException {
        pinMap.clear();
        PinReader pr = new PinReader(pinFile);
        while (pr.next(true)) {
            PeptideFormatter pf = pr.getPep();
            pinMap.put(pr.getScanNum() + "|" + pr.getRank(), pf.getBaseCharge());
        }
    }

    //baseCharge to unicode bin
    public void mapPepxmlToPin() {
        if (pinMap.size() == 0) { //empty
            printError("must load pin before mapping pepxml to pin");
            System.exit(1);
        }

        //get rank
        String rank = "1";
        if (pepxmlName.contains("_rank")) {
            String[] nameSplit = pepxmlName.split("_rank");
            String rankInfo = nameSplit[nameSplit.length - 1];
            rank = rankInfo.split(".pepXML")[0];
        }

        for (int i = 0; i < divisions; i++) {
            for (SpectrumQuery sq : sqToBin[i]) {
                String[] spectrumSplit = sq.getSpectrum().split("\\.");
                String scanNum = spectrumSplit[spectrumSplit.length - 2];

                String baseCharge = pinMap.get(scanNum + "|" + rank);
                baseChargetoUnicodeBin.put(baseCharge, i);
            }
        }
    }

    //alphapeptdeep
    public void inputFileAssignBins(String inputFile) throws IOException {
        //load file
        BufferedReader br = new BufferedReader(new FileReader(inputFile));
        String line = br.readLine();
        String[] header = line.split(",");

        //set column indices
        int chargeIdx = ArrayUtils.indexOf(header, "charge");
        int baseIdx = ArrayUtils.indexOf(header, "base");

        //open {divisions} files and write rows to appropriate files
        BufferedWriter[] writerArray = new BufferedWriter[divisions];
        for (int i = 0; i < writerArray.length; i++) {
            String outFile = inputFile.replace(".csv", i + ".csv");
            writerArray[i] = new BufferedWriter(new FileWriter(outFile));
            writerArray[i].write(line + "\n");
        }

        while ((line = br.readLine()) != null) {
            String[] lineSplit = line.split(",");
            int bin = baseChargetoUnicodeBin.get(lineSplit[baseIdx] + "|" + lineSplit[chargeIdx]);
            writerArray[bin].write(line + "\n");
        }

        for (BufferedWriter bufferedWriter : writerArray) {
            bufferedWriter.close();
        }
    }

    public void dividePinPepxml(File[] pins, String predInputFile) throws IOException, JAXBException, FileParsingException {
        HashMap<String, ArrayList<String>> pinToPepxml = new HashMap<>();
        //collect folders for all pin files
        HashSet<String> folders = new HashSet<>();
        for (File pin : pins) {
            folders.add(pin.getParent());
        }
        for (String folder : folders) {
            File folderFile = new File(folder);
            File[] pepxmls = folderFile.listFiles((dir, name) -> name.endsWith(".pepXML"));
            for (File f : pepxmls) {
                String name = f.getCanonicalPath();
                if (name.contains("divideMSBooster")) {
                    continue;
                }
                if (name.contains("_rank")) {
                    String[] nameSplit = name.split("_rank");
                    if (nameSplit.length > 2) {
                        StringBuilder sb = new StringBuilder();
                        for (int i = 0; i < nameSplit.length - 1; i++) {
                            sb.append(nameSplit[i]);
                            if (i < nameSplit.length - 2) {
                                sb.append("_rank");
                            }
                        }
                        name = sb.toString();
                    } else {
                        name = nameSplit[0];
                    }
                } else {
                    name = name.split(".pepXML")[0];
                }
                name += ".pin";

                ArrayList<String> a;
                if (!pinToPepxml.containsKey(name)) {
                    a = new ArrayList<>();
                } else {
                    a = pinToPepxml.get(name);
                }
                a.add(f.getCanonicalPath());
                pinToPepxml.put(name, a);
            }
        }

        for (File pin : pins) {
            loadPin(pin.getCanonicalPath());

            for (String pepxml : pinToPepxml.get(pin.getCanonicalPath())) {
                printInfo("Processing " + pepxml);
                process(pepxml);
                mapPepxmlToPin();
            }
        }

        printInfo("Creating divided AlphaPeptDeep input files");
        inputFileAssignBins(predInputFile);
    }
}
