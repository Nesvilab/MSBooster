import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsPipelineAnalysis;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsRunSummary;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SearchSummary;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SpectrumQuery;

import javax.xml.bind.JAXBContext;
import javax.xml.bind.JAXBException;
import javax.xml.bind.Marshaller;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.lang.reflect.Method;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.Stream;

public class pepXMLModifier {

    final pepXMLReader xmlReader;
    final String[] xmlPeptides;

    public pepXMLModifier(String pepXMLFile) {
        //create new pepXMLReader object
        xmlReader = new pepXMLReader(pepXMLFile);

        //get peptide names compatible with predicted HashMaps
        xmlPeptides = xmlReader.getXMLpeptides();
    }

    public MsmsPipelineAnalysis replaceExpectScore(mzMLReader mzMLscans, mgfFileReader predictedSpectra,
                                                   String similarityMeasure)
            throws IOException, FileParsingException, NoSuchMethodException,
            InvocationTargetException, IllegalAccessException {

        //start indexing xmlPeptides
        int xmlPeptidesIndex = -1;

        //only use if need weighted similarity
        boolean useWeights = false;
        double binwidth = 1;
        double[] mzFreqs = new double[0];
        if (similarityMeasure.substring(0, 8).equals("weighted")) {
            mzFreqs = mzMLscans.getMzFreq(binwidth, 1);
            useWeights = true;
        }

        //get predictions for unmodified spectra. Skip if modified for now
        HashMap<String, double[]> predictedMZs = predictedSpectra.getMzDict();
        HashMap<String, double[]> predictedIntensities = predictedSpectra.getIntensityDict();

        //the real stuff
        //make new pepxml
        MsmsRunSummary oldRunSummary = xmlReader.analysis.getMsmsRunSummary().get(0);

        MsmsPipelineAnalysis newMsmsPipelineAnalysis = new MsmsPipelineAnalysis();
        List<MsmsRunSummary> newMsmsRunSummaryList = newMsmsPipelineAnalysis.getMsmsRunSummary();
        MsmsRunSummary newMsmsRunSummary = new MsmsRunSummary();
        newMsmsRunSummaryList.add(newMsmsRunSummary);

        //edit MSMSRunSummary to match old one
        newMsmsRunSummary.setBaseName(oldRunSummary.getBaseName());
        newMsmsRunSummary.setRawData(oldRunSummary.getRawData());
        newMsmsRunSummary.setRawDataType(oldRunSummary.getRawDataType());
        newMsmsRunSummary.setSampleEnzyme(oldRunSummary.getSampleEnzyme());
        List<SearchSummary> searchSummaries = newMsmsRunSummary.getSearchSummary();
        searchSummaries.addAll(oldRunSummary.getSearchSummary());

        //modify score
        List<SpectrumQuery> newSpectrumQueryList = newMsmsRunSummary.getSpectrumQuery();
        List<SpectrumQuery> spectrumQueries = oldRunSummary.getSpectrumQuery();
        for (SpectrumQuery sq : spectrumQueries) {
            //get predicted info
            xmlPeptidesIndex += 1;
            String pep = xmlPeptides[xmlPeptidesIndex];
            double[] predMZs = predictedMZs.get(pep);
            double[] predIntensities = predictedIntensities.get(pep);

            //ignore U amino acids
            if (predMZs == null) {
                continue;
            }

            //get experimental info
            //get scan num and query mzml
            int scanNum = (int) sq.getStartScan();
            double[] expMZs = mzMLscans.getMZ(scanNum);
            double[] expIntensities = mzMLscans.getIntensity(scanNum);

            //spectral similarity calculation
            //multiply by negative one, since for log(escore), smaller is better
            spectrumComparison specAngle = new spectrumComparison(expMZs, expIntensities,
                    predMZs, predIntensities, 20);

            //check that similarityMeasure is valid

            double sim = 0;

            //for not getting similarities of zero
            double correctionFactor = 0.1;
            if (similarityMeasure.equals("euclideanDistance") ||
                    similarityMeasure.equals("weightedEuclideanDistance")) {
                correctionFactor += Math.sqrt(2) - 1;
            }

            //only if need weights
            if (useWeights) {
                Method method = specAngle.getClass().getMethod(similarityMeasure, double[].class);
                double[] weights = specAngle.getWeights(mzFreqs, binwidth);
                sim = (double) method.invoke(specAngle, weights);
            } else {
                Method method = specAngle.getClass().getMethod(similarityMeasure);
                sim = (double) method.invoke(specAngle);
            }

            sim = 1 / (sim + correctionFactor);

            //can we make assumption there's only 1 element per list?
            //set new substitute for expectation value
            sq.getSearchResult().get(0).getSearchHit().get(0).getSearchScore().get(2).setValueStr(String.valueOf(sim));
            newSpectrumQueryList.add(sq);
        }
        return newMsmsPipelineAnalysis;
    }

    public MsmsPipelineAnalysis hyperscoreXML()
            throws IOException, FileParsingException, NoSuchMethodException,
            InvocationTargetException, IllegalAccessException {

        //the real stuff
        //make new pepxml
        MsmsRunSummary oldRunSummary = xmlReader.analysis.getMsmsRunSummary().get(0);

        MsmsPipelineAnalysis newMsmsPipelineAnalysis = new MsmsPipelineAnalysis();
        List<MsmsRunSummary> newMsmsRunSummaryList = newMsmsPipelineAnalysis.getMsmsRunSummary();
        MsmsRunSummary newMsmsRunSummary = new MsmsRunSummary();
        newMsmsRunSummaryList.add(newMsmsRunSummary);

        //edit MSMSRunSummary to match old one
        newMsmsRunSummary.setBaseName(oldRunSummary.getBaseName());
        newMsmsRunSummary.setRawData(oldRunSummary.getRawData());
        newMsmsRunSummary.setRawDataType(oldRunSummary.getRawDataType());
        newMsmsRunSummary.setSampleEnzyme(oldRunSummary.getSampleEnzyme());
        List<SearchSummary> searchSummaries = newMsmsRunSummary.getSearchSummary();
        searchSummaries.addAll(oldRunSummary.getSearchSummary());

        //modify score
        List<SpectrumQuery> newSpectrumQueryList = newMsmsRunSummary.getSpectrumQuery();
        List<SpectrumQuery> spectrumQueries = oldRunSummary.getSpectrumQuery();
        for (SpectrumQuery sq : spectrumQueries) {
            //can we make assumption there's only 1 element per list?
            //set new substitute for expectation value
            String hyperscoreStr = sq.getSearchResult().get(0).getSearchHit().get(0).getSearchScore().get(0).getValueStr();
            double hyperscore = Double.parseDouble(hyperscoreStr);
            hyperscore = hyperscore = 10 / (hyperscore + 0.1);

            sq.getSearchResult().get(0).getSearchHit().get(0).getSearchScore().get(2).setValueStr(String.valueOf(hyperscore));
            newSpectrumQueryList.add(sq);
        }
        return newMsmsPipelineAnalysis;
    }

    public void writepepXML(MsmsPipelineAnalysis newMsmsPipelineAnalysis, String outFile) throws JAXBException {
        //write xml
        JAXBContext jaxbContext = JAXBContext.newInstance(MsmsPipelineAnalysis.class);
        Marshaller marshaller = jaxbContext.createMarshaller();
        marshaller.setProperty(Marshaller.JAXB_FORMATTED_OUTPUT, Boolean.TRUE);

        Path path = Paths.get(outFile);
        java.io.File file = path.toFile();
        marshaller.marshal(newMsmsPipelineAnalysis, file);
    }

    public static void main(String[] args) throws IOException, FileParsingException, NoSuchMethodException, IllegalAccessException, InvocationTargetException, JAXBException {
        for (int i = 1; i < 7; i++) {
            String narrowNum = String.valueOf(i);

            Path start = Paths.get("C:/Users/kevin/OneDriveUmich/proteomics/pepxml/");

            String[] windowFiles = new String[4];
            try (Stream<Path> stream = Files.walk(start)) {
                List<String> collect = stream
                        .map(String::valueOf)
                        .sorted()
                        .collect(Collectors.toList());

                int filesIndex = 0;
                for (String s : collect) {
                    if (s.contains("narrow_" + narrowNum)) {
                        windowFiles[filesIndex] = s;
                        filesIndex += 1;
                    }
                }
            }

            //modify files
            mzMLReader mzmlScans = new mzMLReader("C:/Users/kevin/OneDriveUmich/proteomics/mzml/" +
                    "23aug2017_hela_serum_timecourse_4mz_narrow_" + narrowNum + ".mzML");
            mgfFileReader preds = new mgfFileReader("preds/"); //stays constant for all mzml files
            String similarityMeasure = "weightedEuclideanDistance";

            for (String file : windowFiles) {
                pepXMLModifier xmlReader = new pepXMLModifier(file);
                MsmsPipelineAnalysis modifiedXML = xmlReader.replaceExpectScore(mzmlScans, preds, similarityMeasure);
                //MsmsPipelineAnalysis modifiedXML = xmlReader.hyperscoreXML();

                String file1 = file.split("\\\\")[file.split("\\\\").length - 1];
                String[] underscoreSplits = file1.split("_");
                String rankDir = underscoreSplits[underscoreSplits.length - 1].split("\\.")[0];
                String fileName = String.join("_",
                        Arrays.copyOfRange(underscoreSplits, 0, underscoreSplits.length - 1)) + ".pepXML";
                String newDirectoryName = "C:/Users/kevin/Downloads/proteomics/" + similarityMeasure + "/" + rankDir + "/";
                File newDirectory = new File(newDirectoryName);
                if (! newDirectory.exists()){
                    newDirectory.mkdirs();
                }
                String outFile = newDirectoryName + fileName;
                System.out.println("writing to " + outFile);
                xmlReader.writepepXML(modifiedXML, outFile);
            }
        }
    }
}