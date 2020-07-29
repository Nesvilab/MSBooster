import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.*;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;

public class pepXMLReader {

    String pathname;
    MsmsPipelineAnalysis analysis;
    int[] scanNumbers;
    String[] XMLpeptides;
    int[] targetOrDecoy;

    public pepXMLReader(String path) {
        pathname = path;
        analysis = this.getAnalysis(pathname);
    }

    public MsmsPipelineAnalysis getAnalysis(String path) {
        Path p = Paths.get(path);
        MsmsPipelineAnalysis a = new MsmsPipelineAnalysis();
        try {
            a = PepXmlParser.parse(p);
        } catch (FileParsingException e) {
            e.printStackTrace();
        }
        return a;
    }

    public void createPDeepFile(String outfile) {

        // iterate over the parsed search results
        List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
        for (MsmsRunSummary runSummary : runSummaries) {

            List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();

            //initialize array of peptide hit info
            String[] allHits = new String[spectrumQueries.size()];
            int allHitsIndex = -1;
            for (SpectrumQuery sq : spectrumQueries) {
                allHitsIndex++;

                //get charge and initialize other variables
                String pep = "";
                String modifications = "";
                String charge = String.valueOf(sq.getAssumedCharge());

                List<SearchResult> results = sq.getSearchResult();

                for (SearchResult result : results) {
                    List<SearchHit> hits = result.getSearchHit();

                    for (SearchHit hit : hits) {
                        pep = hit.getPeptide();
                        modifications = "";

                        ModificationInfo modinfo = hit.getModificationInfo();
                        if (modinfo != null) {
                            //n term acetylation
                            if (modinfo.getModNtermMass() != null) {
                                modifications += "0,Acetyl[AnyN-term];";
                            }

                            //modifications besides n term acetylation
                            List<ModAminoacidMass> mods = modinfo.getModAminoacidMass();
                            int modSize = mods.size();
                            if (modSize != 0) {
                                int[] pos = new int[modSize];
                                double[] modMass = new double[modSize];

                                int i = 0;
                                for (ModAminoacidMass mod : mods) {
                                    pos[i] = mod.getPosition();
                                    modMass[i] = mod.getMass();
                                    i ++;
                                }

                                //format new modifications
                                modFormatter mf = new modFormatter(pos, modMass);

                                //add formatted new modifications onto modifications string
                                modifications += mf.format();
                            }
                        }
                    }
                }
                String hitToAdd = pep + "\t" + modifications + "\t" + charge;
                allHits[allHitsIndex] = hitToAdd;
            }

            //remove duplicates from allHits
            //can reduce number of hits to a third
            HashSet<String> hSetHits = new HashSet<>();
            Collections.addAll(hSetHits, allHits);

            //write to file
            try {
                FileWriter myWriter = new FileWriter(outfile);
                myWriter.write("peptide" + "\t" + "modification" + "\t" + "charge\n");

                for (String hSetHit : hSetHits) {
                    myWriter.write(hSetHit + "\n");
                }

                myWriter.close();
                System.out.println("Successfully wrote to the file.");
            } catch (IOException e) {
                System.out.println("An error occurred.");
                e.printStackTrace();
            }
        }
    }

    //method to return iterable of scan numbers
    public int[] getScanNumbers() {
        if (this.scanNumbers == null) {
            int[] scanNumbers = new int[0];
            List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
            for (MsmsRunSummary runSummary : runSummaries) {
                List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();

                int querySize = spectrumQueries.size();
                scanNumbers = new int[querySize];
                for (int i = 0; i < querySize; i++) {
                    scanNumbers[i] = (int) spectrumQueries.get(i).getStartScan();
                }
            }
            this.scanNumbers = scanNumbers;
            return scanNumbers;
        } else {
            return this.scanNumbers;
        }
    }

    public String[] getXMLpeptides() {
        if (XMLpeptides == null) {

            String[] peptides = new String[0];

            List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
            for (MsmsRunSummary runSummary : runSummaries) {
                List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();
                //initialize array of peptide hit info
                peptides = new String[spectrumQueries.size()];
                int index = -1;
                for (SpectrumQuery sq : spectrumQueries) {
                    index++;

                    //get charge and initialize other variables
                    String pep = "";
                    String modifications = "";
                    String charge = String.valueOf(sq.getAssumedCharge());

                    List<SearchResult> results = sq.getSearchResult();

                    for (SearchResult result : results) {
                        List<SearchHit> hits = result.getSearchHit();

                        for (SearchHit hit : hits) {
                            pep = hit.getPeptide();
                            modifications = "";

                            ModificationInfo modinfo = hit.getModificationInfo();
                            if (modinfo != null) {
                                //n term acetylation
                                if (modinfo.getModNtermMass() != null) {
                                    modifications += "0,Acetyl[AnyN-term];";
                                }

                                //modifications besides n term acetylation
                                List<ModAminoacidMass> mods = modinfo.getModAminoacidMass();
                                int modSize = mods.size();
                                if (modSize != 0) {
                                    int[] pos = new int[modSize];
                                    double[] modMass = new double[modSize];

                                    int i = 0;
                                    for (ModAminoacidMass mod : mods) {
                                        pos[i] = mod.getPosition();
                                        modMass[i] = mod.getMass();
                                        i++;
                                    }

                                    //format new modifications
                                    modFormatter mf = new modFormatter(pos, modMass);

                                    //add formatted new modifications onto modifications string
                                    modifications += mf.format();
                                }
                            }
                        }
                    }

                    String XMLpeptide = pep + "|" + modifications + "|" + charge;
                    peptides[index] = XMLpeptide;
                }
            }
            this.XMLpeptides = peptides;
            return peptides;
        } else {
            return this.XMLpeptides;
        }
    }

    //method to return iterable of scan numbers
    public int[] getTargetOrDecoy() {
        if (this.targetOrDecoy == null) {
            int[] td = new int[0];

            List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
            for (MsmsRunSummary runSummary : runSummaries) {
                List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();
                int i = 0;
                td = new int[spectrumQueries.size()];

                for (SpectrumQuery sq : spectrumQueries) {
                    List<SearchResult> results = sq.getSearchResult();
                    for (SearchResult result : results) {
                        List<SearchHit> hits = result.getSearchHit();
                        for (SearchHit hit : hits) {
                            String prefix = hit.getProtein().substring(0, 3);
                            if (prefix.equals("rev")) {
                                td[i] = 0; //decoy
                            } else {
                                td[i] = 1; //target
                            }
                            i++;
                        }
                    }
                }
            }
            this.targetOrDecoy = td;
            return td;
        } else {
            return this.targetOrDecoy;
        }
    }

    public static void main(String[] args) {
        pepXMLReader x = new pepXMLReader("23aug2017_hela_serum_timecourse_4mz_narrow_1_rank1.pepXML");
        List<MsmsRunSummary> runSummaries = x.analysis.getMsmsRunSummary();
        System.out.println(runSummaries.size());
    }
}
