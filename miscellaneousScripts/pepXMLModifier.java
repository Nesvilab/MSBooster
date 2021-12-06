package Features;

import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsPipelineAnalysis;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.MsmsRunSummary;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SearchSummary;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.SpectrumQuery;

import javax.xml.bind.JAXBException;
import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
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
        double[] mzFreqs = new double[0];
        if (similarityMeasure.substring(0, 6).equals("weight")) {
            mzFreqs = mzMLscans.getMzFreq();
            useWeights = true;
        }

        //get predictions for unmodified spectra. Skip if modified for now
        HashMap<String, float[]> predictedMZs = predictedSpectra.getMzDict();
        HashMap<String, float[]> predictedIntensities = predictedSpectra.getIntensityDict();

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
            float[] predMZs = predictedMZs.get(pep);
            float[] predIntensities = predictedIntensities.get(pep);

            //ignore U amino acids
            if (predMZs == null) {
                continue;
            }

            //get experimental info
            //get scan num and query mzml
            int scanNum = (int) sq.getStartScan();
            float[] expMZs = mzMLscans.getMZ(scanNum);
            float[] expIntensities = mzMLscans.getIntensity(scanNum);

            //spectral similarity calculation
            //multiply by negative one, since for log(escore), smaller is better
            spectrumComparison specAngle = new spectrumComparison(expMZs, expIntensities,
                    predMZs, predIntensities, Constants.useTopFragments);

            //check that similarityMeasure is valid

            //for not getting similarities of zero
            double correctionFactor = 0.1;
            if (similarityMeasure.equals("euclideanDistance") ||
                    similarityMeasure.equals("weightedEuclideanDistance")) {
                correctionFactor += Math.sqrt(2) - 1;
            }

            double sim = specAngle.getSimilarity(similarityMeasure, mzFreqs);

            //if getting expectation score, calculate here
            //create scannumobjects for mzmlreader
            //then using scannum, retrieve mzmlscannumber instance from scannumberobjects, and get its windowstart field
            //get other scannums, and
            //get window start from

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

    public static void main(String[] args) throws IOException, FileParsingException, NoSuchMethodException,
            IllegalAccessException, InvocationTargetException, JAXBException {
        for (int i = 1; i < 7; i++) {
            String narrowNum = String.valueOf(i);

            Path start = Paths.get("C:/Users/kevin/Downloads/proteomics/greedyDoubleCounting/OGrank");

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
            mgfFileReader preds = new mgfFileReader("C:/Users/kevin/Downloads/proteomics/" +
                    "pDeep3preds_new.mgf"); //stays constant for all mzml files
            String similarityMeasure = "brayCurtis";

            for (int j = 0; j < windowFiles.length; j++) {
                String file = windowFiles[j];
                pepXMLModifier xmlReader = new pepXMLModifier(file);
                MsmsPipelineAnalysis modifiedXML = xmlReader.replaceExpectScore(mzmlScans, preds, similarityMeasure);
                //MsmsPipelineAnalysis modifiedXML = xmlReader.hyperscoreXML();

                String file1 = file.split("\\\\")[file.split("\\\\").length - 1];
                String newDirectoryName = "C:/Users/kevin/Downloads/proteomics/greedyDoubleCounting/"
                        + similarityMeasure + "/rank" + String.valueOf(j + 1) + "/";
                File newDirectory = new File(newDirectoryName);
                if (! newDirectory.exists()){
                    newDirectory.mkdirs();
                }
                String outFile = newDirectoryName + file1;
                System.out.println("writing to " + outFile);
                pepXMLWriter.writepepXML(modifiedXML, outFile);
            }
        }
    }
}