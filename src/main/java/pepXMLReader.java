import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.*;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class pepXMLReader {

    String pathname;
    MsmsPipelineAnalysis analysis;
    int[] scanNumbers;
    String[] XMLpeptides;
    int[] targetOrDecoy;
    String[] allHitsPDeep;
    String[] eScores;

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

    public void createPDeepList() {

        // iterate over the parsed search results
        List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
        for (MsmsRunSummary runSummary : runSummaries) {

            List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();

            //initialize array of peptide hit info
            allHitsPDeep = new String[spectrumQueries.size()];
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
                allHitsPDeep[allHitsIndex] = hitToAdd;
            }
        }
    }

    public void createPDeepListNoMods() {

        // iterate over the parsed search results
        List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
        for (MsmsRunSummary runSummary : runSummaries) {

            List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();

            //initialize array of peptide hit info
            allHitsPDeep = new String[spectrumQueries.size()];
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

                            //modifications besides n term acetylation
                            List<ModAminoacidMass> mods = modinfo.getModAminoacidMass();
                            for (int i = 0; i < mods.size(); i++) {
                                ModAminoacidMass mod = mods.get(i);
                                if (mod.getMass() == 147.0354) {
                                    mods.remove(i);
                                }
                            }

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
                allHitsPDeep[allHitsIndex] = hitToAdd;
            }
        }
    }

    public void createPDeep3List() {

        // iterate over the parsed search results
        List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
        for (MsmsRunSummary runSummary : runSummaries) {

            List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();

            //initialize array of peptide hit info
            allHitsPDeep = new String[spectrumQueries.size()];
            int allHitsIndex = -1;
            for (SpectrumQuery sq : spectrumQueries) {
                allHitsIndex++;

                //get charge and initialize other variables
                String pep = "";
                String modifications = "";
                String charge = String.valueOf(sq.getAssumedCharge());
                //int scanNum = (int) sq.getStartScan();
                //String scanSumStr = String.valueOf(scanNum);

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
                //I don't think raw name or scanNum is used for predictions?
                String hitToAdd = "." + "\t" + "." + "\t" + pep + "\t" + modifications + "\t" + charge;
                allHitsPDeep[allHitsIndex] = hitToAdd;
            }
        }
    }

    public void createPrositList(int NCE) {

        // iterate over the parsed search results
        List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
        for (MsmsRunSummary runSummary : runSummaries) {

            List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();

            //initialize array of peptide hit info
            allHitsPDeep = new String[spectrumQueries.size()];
            int allHitsIndex = -1;
            for (SpectrumQuery sq : spectrumQueries) {
                allHitsIndex++;

                //get charge and initialize other variables
                String pep = "";
                String charge = String.valueOf(sq.getAssumedCharge());
                boolean skip = false;

                List<SearchResult> results = sq.getSearchResult();

                for (SearchResult result : results) {
                    List<SearchHit> hits = result.getSearchHit();

                    for (SearchHit hit : hits) {
                        pep = hit.getPeptide();
                        if (pep.length() > 27) {
                            skip = true;
                            continue;
                        }

                        ModificationInfo modinfo = hit.getModificationInfo();
                        if (modinfo != null) {
                            //n term acetylation
                            if (modinfo.getModNtermMass() != null) {
                                //Prosit does not support n-term acetylation
                                skip = true;
                            } else { //safe to continue
                                //modifications besides n term acetylation
                                List<ModAminoacidMass> mods = modinfo.getModAminoacidMass();
                                int modSize = mods.size();
                                if (modSize != 0) {
                                    for (ModAminoacidMass mod : mods) {
                                        if (mod.getMass() == 147.0354) { //just check for oxidized Met
                                            pep = pep.substring(0, mod.getPosition()) + 'm' +
                                                    pep.substring(mod.getPosition());
                                        }
                                    }
                                    pep = pep.replace("m", "M(ox)");
                                }
                            }
                        }
                    }
                }
                if (! skip) {
                    String hitToAdd = pep + "," + NCE + "," + charge; //HeLa run with 27 NCE
                    allHitsPDeep[allHitsIndex] = hitToAdd;
                }
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

    public String[] getEScore() {
        if (this.eScores == null) {
            String[] e = new String[0];

            List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
            for (MsmsRunSummary runSummary : runSummaries) {
                List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();
                int i = 0;
                e = new String[spectrumQueries.size()];

                for (SpectrumQuery sq : spectrumQueries) {
                    List<SearchResult> results = sq.getSearchResult();
                    for (SearchResult result : results) {
                        List<SearchHit> hits = result.getSearchHit();
                        for (SearchHit hit : hits) {
                            e[i] = hit.getSearchScore().get(2).getValueStr();
                            i++;
                        }
                    }
                }
            }
            this.eScores = e;
            return e;
        } else {
            return this.eScores;
        }
    }

    public static void main(String[] args) {
        pepXMLReader x = new pepXMLReader("C:/Users/kevin/Downloads/proteomics/expect/rank1/" +
                "23aug2017_hela_serum_timecourse_4mz_narrow_1_rank1.pepXML");
        System.out.println(x.getScanNumbers().length);

    }
}


