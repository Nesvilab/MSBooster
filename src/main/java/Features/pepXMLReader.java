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

package Features;

import umich.ms.fileio.exceptions.FileParsingException;
import umich.ms.fileio.filetypes.pepxml.PepXmlParser;
import umich.ms.fileio.filetypes.pepxml.jaxb.standard.*;

import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.List;

public class pepXMLReader {

    String pathname;
    MsmsPipelineAnalysis analysis;
    String baseName;
    private int[] scanNumbers; //private means I have to use get{something} method
    private String[] XMLpeptides;
    private int[] targetOrDecoy;
    private String[] eScores;
    private double[] RTs;
    private double[] massDiffs;
    private int[] NTTs;
    private int[] NMCs;
    private int[] isotopeErrs;

    public pepXMLReader(String path) {
        pathname = path;
        analysis = this.getAnalysis(pathname);
        String[] baseNameArray = analysis.getMsmsRunSummary().get(0).getBaseName().split("\\\\");
        baseName = baseNameArray[baseNameArray.length - 1];
        //baseName = analysis.getMsmsRunSummary().get(0).getBaseName();
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

//    public void createPDeepList() {
//
//        // iterate over the parsed search results
//        List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
//        for (MsmsRunSummary runSummary : runSummaries) {
//
//            List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();
//
//            //initialize array of peptide hit info
//            allHitsPDeep = new String[spectrumQueries.size()];
//            int allHitsIndex = -1;
//            for (SpectrumQuery sq : spectrumQueries) {
//                allHitsIndex++;
//
//                //get charge and initialize other variables
//                String pep = "";
//                String modifications = "";
//                String charge = String.valueOf(sq.getAssumedCharge());
//
//                List<SearchResult> results = sq.getSearchResult();
//
//                for (SearchResult result : results) {
//                    List<SearchHit> hits = result.getSearchHit();
//
//                    for (SearchHit hit : hits) {
//                        pep = hit.getPeptide();
//                        modifications = "";
//
//                        ModificationInfo modinfo = hit.getModificationInfo();
//                        if (modinfo != null) {
//                            //n term acetylation
//                            if (modinfo.getModNtermMass() != null) {
//                                modifications += "0,Acetyl[AnyN-term];";
//                            }
//
//                            //modifications besides n term acetylation
//                            List<ModAminoacidMass> mods = modinfo.getModAminoacidMass();
//                            int modSize = mods.size();
//                            if (modSize != 0) {
//                                int[] pos = new int[modSize];
//                                double[] modMass = new double[modSize];
//
//                                int i = 0;
//                                for (ModAminoacidMass mod : mods) {
//                                    pos[i] = mod.getPosition();
//                                    modMass[i] = mod.getMass();
//                                    i ++;
//                                }
//
//                                //format new modifications
//                                modFormatter mf = new modFormatter(pos, modMass);
//
//                                //add formatted new modifications onto modifications string
//                                modifications += mf.format();
//                            }
//                        }
//                    }
//                }
//                String hitToAdd = pep + "\t" + modifications + "\t" + charge;
//                allHitsPDeep[allHitsIndex] = hitToAdd;
//            }
//        }
//    }

//    public void createPDeepListNoMods() {
//
//        // iterate over the parsed search results
//        List<MsmsRunSummary> runSummaries = analysis.getMsmsRunSummary();
//        for (MsmsRunSummary runSummary : runSummaries) {
//
//            List<SpectrumQuery> spectrumQueries = runSummary.getSpectrumQuery();
//
//            //initialize array of peptide hit info
//            allHitsPDeep = new String[spectrumQueries.size()];
//            int allHitsIndex = -1;
//            for (SpectrumQuery sq : spectrumQueries) {
//                allHitsIndex++;
//
//                //get charge and initialize other variables
//                String pep = "";
//                String modifications = "";
//                String charge = String.valueOf(sq.getAssumedCharge());
//
//                List<SearchResult> results = sq.getSearchResult();
//
//                for (SearchResult result : results) {
//                    List<SearchHit> hits = result.getSearchHit();
//
//                    for (SearchHit hit : hits) {
//                        pep = hit.getPeptide();
//                        modifications = "";
//
//                        ModificationInfo modinfo = hit.getModificationInfo();
//                        if (modinfo != null) {
//
//                            //modifications besides n term acetylation
//                            List<ModAminoacidMass> mods = modinfo.getModAminoacidMass();
//                            for (int i = 0; i < mods.size(); i++) {
//                                ModAminoacidMass mod = mods.get(i);
//                                if (mod.getMass() == 147.0354) {
//                                    mods.remove(i);
//                                }
//                            }
//
//                            int modSize = mods.size();
//                            if (modSize != 0) {
//                                int[] pos = new int[modSize];
//                                double[] modMass = new double[modSize];
//
//                                int i = 0;
//                                for (ModAminoacidMass mod : mods) {
//                                    pos[i] = mod.getPosition();
//                                    modMass[i] = mod.getMass();
//                                    i ++;
//                                }
//
//                                //format new modifications
//                                modFormatter mf = new modFormatter(pos, modMass);
//
//                                //add formatted new modifications onto modifications string
//                                modifications += mf.format();
//                            }
//                        }
//                    }
//                }
//                String hitToAdd = pep + "\t" + modifications + "\t" + charge;
//                allHitsPDeep[allHitsIndex] = hitToAdd;
//            }
//        }
//    }

    public String[] createPDeep3List() {

        List<SpectrumQuery> spectrumQueries = getSpectrumQueries();

        //initialize array of peptide hit info
        String[] allHitsPDeep = new String[spectrumQueries.size()];
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
        return allHitsPDeep;
    }

    public String[] createPrositList(int NCE) {
        List<SpectrumQuery> spectrumQueries = getSpectrumQueries();
        String[] array = new String[spectrumQueries.size()];

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
                    if (pep.length() > 30) {
                        skip = true;
                        break;
                    }
                    if (pep.contains("O") || pep.contains("U") || pep.contains("J")) {
                        skip = true;
                        break;
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
                                        pep = pep.substring(0, mod.getPosition() - 1) + 'm' +
                                                pep.substring(mod.getPosition());
                                    }
                                }
                                pep = pep.replace("m", "M(ox)");
                            }
                        }
                    }
                }
            }
            if (!skip) {
                String hitToAdd = pep + "," + NCE + "," + charge; //HeLa run with 27 NCE
                array[allHitsIndex] = hitToAdd;
            }
        }
        return array;
    }

    public String[] createDeepMSPeptideList(){
        List<SpectrumQuery> spectrumQueries = getSpectrumQueries();
        String[] array = new String[spectrumQueries.size()];
        int index = 0;
        for (SpectrumQuery sq : spectrumQueries) {
            array[index] = sq.getSearchResult().get(0).getSearchHit().get(0).getPeptide();
            index++;
        }
        return array;
    }

//    public String[] createDiannList() {
//        List<SpectrumQuery> spectrumQueries = getSpectrumQueries();
//        String[] array = new String[spectrumQueries.size()];
//        int index = 0;
//        for (SpectrumQuery sq : spectrumQueries) {
//            String charge = String.valueOf(sq.getAssumedCharge());
//            SearchHit hit = sq.getSearchResult().get(0).getSearchHit().get(0);
//            StringBuilder pep = new StringBuilder(hit.getPeptide());
//
//            ModificationInfo modinfo = hit.getModificationInfo();
//            TreeMap<Integer, Integer> modMap = new TreeMap<>(); //sorted for future use
//            if (modinfo != null) {
//                //n term acetylation
//                if (modinfo.getModNtermMass() != null) {
//                    modMap.put(0, Constants.modNtermToUnimod.get(modinfo.getModNtermMass()));
//                }
//
//                //modifications besides n term acetylation
//                List<ModAminoacidMass> mods = modinfo.getModAminoacidMass();
//                for (ModAminoacidMass mod : mods) {
//                    //match modmass to unimod
//                    //pep.insert(mod.getPosition(), "[unimod:" + constants.modAAmassToUnimod.get(mod.getMass()) + "]");
//                    modMap.put(mod.getPosition(), Constants.modAAmassToUnimod.get(mod.getMass()));
//                }
//
//                //add mods to peptide
//                int strLen = 0;
//                for (Map.Entry<Integer, Integer> entry : modMap.entrySet()) {
//                    String mod = "[unimod:" + entry.getValue() + "]";
//                    pep.insert(strLen + entry.getKey(), mod);
//                    strLen += mod.length();
//                }
//
//            }
//
//            //set peptide
//            array[index] = pep.toString() + "\t" + charge;
//            index++;
//        }
//        return array;
//    }

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

//    public static int[] getScanNumbers(List<SpectrumQuery> spectrumQueries) {
//        int querySize = spectrumQueries.size();
//
//        int[] scanNumbers = new int[querySize];
//        for (int i = 0; i < querySize; i++) {
//            scanNumbers[i] = (int) spectrumQueries.get(i).getStartScan();
//        }
//
//        return scanNumbers;
//    }

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

    public double[] getRT() {
        if (this.RTs == null) {
            List<SpectrumQuery> sqs = analysis.getMsmsRunSummary().get(0).getSpectrumQuery();
            double[] rt = new double[sqs.size()];
            int i = 0;
            for (SpectrumQuery sq : sqs) {
                rt[i] = sq.getRetentionTimeSec();
                i++;
            }

            this.RTs = rt;
            return rt;
        } else {
            return this.RTs;
        }
    }

    public double[] getMassDiff() {
        if (this.massDiffs == null) {
            List<SpectrumQuery> sqs = analysis.getMsmsRunSummary().get(0).getSpectrumQuery();
            double[] md = new double[sqs.size()];
            int i = 0;
            for (SpectrumQuery sq : sqs) {
                md[i] = sq.getSearchResult().get(0).getSearchHit().get(0).getMassdiff();
                i++;
            }

            this.massDiffs = md;
            return md;
        } else {
            return this.massDiffs;
        }
    }

    public int[] getNTT() {
        if (this.NTTs == null) {
            List<SpectrumQuery> sqs = analysis.getMsmsRunSummary().get(0).getSpectrumQuery();
            int[] ntt = new int[sqs.size()];
            int i = 0;
            for (SpectrumQuery sq : sqs) {
                ntt[i] = sq.getSearchResult().get(0).getSearchHit().get(0).getNumTolTerm();
                i++;
            }

            this.NTTs = ntt;
            return ntt;
        } else {
            return this.NTTs;
        }
    }

    public int[] getNMCs() {
        if (this.NMCs == null) {
            List<SpectrumQuery> sqs = analysis.getMsmsRunSummary().get(0).getSpectrumQuery();
            int[] nmc = new int[sqs.size()];
            int i = 0;
            for (SpectrumQuery sq : sqs) {
                nmc[i] = sq.getSearchResult().get(0).getSearchHit().get(0).getNumMissedCleavages();
                i++;
            }

            this.NMCs = nmc;
            return nmc;
        } else {
            return this.NMCs;
        }
    }

    public int[] getIsotopeErrs() {
        if (this.isotopeErrs == null) {
            double[] md = this.getMassDiff();

            int[] ie = new int[md.length];
            for (int i = 0; i < md.length; i++) {
                ie[i] = (int) Math.round(Math.abs(md[i] / 1.0033548378));
            }
            this.isotopeErrs = ie;
            return ie;
        } else {
            return this.isotopeErrs;
        }
    }

    //could be used in other get methods
    public List<SpectrumQuery> getSpectrumQueries() {
        return analysis.getMsmsRunSummary().get(0).getSpectrumQuery();
    }

    public MsmsPipelineAnalysis freshMsmsPipelineAnalysis() {
        MsmsRunSummary oldRunSummary = analysis.getMsmsRunSummary().get(0);

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

        //set spectrum query list. This is where you add new spectrum queries (PSM entries)
        List<SpectrumQuery> newSpectrumQueryList = newMsmsRunSummary.getSpectrumQuery();

        return newMsmsPipelineAnalysis;
    }
}


