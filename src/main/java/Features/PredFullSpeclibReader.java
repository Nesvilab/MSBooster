package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;

public class PredFullSpeclibReader extends mgfFileReader{
    public PredFullSpeclibReader(String file, boolean createScanNumObjects, ExecutorService executorService)
            throws InterruptedException, ExecutionException, FileParsingException, IOException {
        //initially start with mgf reading
        super(file, createScanNumObjects, executorService);

        //also need to add to allPreds the peptides with PTMs not supported by PredFull
        //this requires that mgf file name has same prefix as full input file
        String fullFile = file.substring(0, file.length() - 4) + "_full.tsv";

        //for each peptide, check if base contain PTMs not supported by predfull
        //if so, need to start shifting
        //can do this by enumerating all possible fragment m/z (method in mass calculator)
        BufferedReader TSVReader = new BufferedReader(new FileReader(fullFile));
        String l;

        while ((l = TSVReader.readLine()) != null) {
            //find locations of PTMs
            ArrayList<Integer> starts = new ArrayList<>();
            ArrayList<Integer> ends = new ArrayList<>();
            for (int i = 0; i < l.length(); i++) {
                if (l.charAt(i) == '[') {
                    starts.add(i);
                } else if (l.charAt(i) == ']') {
                    ends.add(i);
                }
            }
            boolean newEntryNeeded = false;
            for (int i = 0; i < starts.size(); i++) {
                Double reportedMass = Double.parseDouble(l.substring(starts.get(i) + 1, ends.get(i)));
                if (!(reportedMass.equals(Constants.carbamidomethylationMass) || reportedMass.equals(Constants.oxidationMass))) {
                    newEntryNeeded = true;
                    break;
                }
            }

            //given that C is assumed to be 57, need to check if no mod on it
            if (! newEntryNeeded) {
                for (int i = 0; i < l.length(); i++) {
                    if (l.charAt(i) == 'C') {
                        if (l.charAt(i + 1) != '[') {
                            newEntryNeeded = true;
                        }
                    }
                }
            }

            if (newEntryNeeded) {
                //start shifting and annotating
                //first, get unmodified peptide prediction from dictionary
                String[] lSplit = l.split("\t");
                PeptideFormatter peptideToSearch = new PeptideFormatter(
                        new PeptideFormatter(lSplit[0], lSplit[1], "base").predfull,
                        lSplit[1], "predfull");

                PredictionEntry pe = allPreds.get(peptideToSearch.baseCharge);
                float[] mzs = pe.mzs;
                float[] intensities = pe.intensities;

                //for each fragment ion, see if it can be annotated
                MassCalculator mc = new MassCalculator(peptideToSearch.base, peptideToSearch.charge);
                mc.possibleFragmentIons();

                //if useBasePeak, this is a quick O(n) process that can limit the number of fragment ions to consider
                double minIntensity = 0f;
                if (Constants.useBasePeak) {
                    for (float i : intensities) {
                        if (i > minIntensity) {
                            minIntensity = i;
                        }
                    }
                    minIntensity = minIntensity * Constants.percentBasePeak / 100f;
                }

                //for each predicted peak, see if it can be annotated
                ArrayList<Float> newMZs = new ArrayList<>();
                ArrayList<Float> newIntensities = new ArrayList<>();
                ArrayList<Float> finalMZs = new ArrayList<>();
                ArrayList<Float> finalIntensities = new ArrayList<>();

                for (int i = 0; i < mzs.length; i++) {
                    if (intensities[i] > minIntensity) {
                        newIntensities.add(intensities[i]);
                        newMZs.add(mzs[i]);
                    }
                }

                //annotate ion
                //any calculated annotated fragment is considered
                String[] annotations = mc.annotateMZs(newMZs);
                //System.out.println(Arrays.toString(annotations));

                //for annotated fragment ions, shift their predicted peak
                //shift only if including the amino acid
                //create new massCalculator object, and change m/z to new one
                //give String in annotation to calcMass and round to closest tenth of m/z
                //adjust in newMZs
                MassCalculator shiftedMC = new MassCalculator(lSplit[0], lSplit[1]);
                int index = 0;
                for (String fragment : annotations) {
                    Set<Float> newMZ = new HashSet<>();

                    //possibly multiple identities, and only accept fragment if all or none are shifted
                    String[] fragments = fragment.split(";");

                    for (String frag : fragments) {
                        //split by + to get change. If no split, assume charge 1
                        String[] chargeSplit = frag.split("\\+");
                        int charge = 1;
                        if (chargeSplit.length == 2) {
                            charge = Integer.parseInt(chargeSplit[1]);
                        }

                        //split by - to get neutral loss. If no split, assume no neutral loss
                        String[] nlSplit = chargeSplit[0].split("-");
                        String nl = "";
                        if (nlSplit.length == 2) {
                            nl = nlSplit[1];
                        }

                        //switch case (MH, immonium, extract letters and numbers)
                        String ionName = nlSplit[0].split(":")[0];
                        switch (ionName) {
                            case "MH":
                                newMZ.add(shiftedMC.calcMass(shiftedMC.modMasses.size() - 2, "y", charge, nl));
                                break;
                            case "immonium":
                                //need to check if mox of modified and unmodified, and exclude if so
                                //don't know if both mod and unmod m/z will be that intensity, safer to exclude
                                Set<Double> mods = new HashSet<>();
                                char aa = frag.charAt(frag.length() - 1);
                                for (int i = 0; i < shiftedMC.peptide.length(); i++) {
                                    if (shiftedMC.peptide.charAt(i) == aa) {
                                        mods.add(shiftedMC.modMasses.get(i + 1));
                                    }
                                }
                                if (mods.size() == 1) {
                                    newMZ.add((float) (mods.iterator().next() + shiftedMC.AAmap.get(aa)));
                                }
                                break;
                            case "unknown":
                                newMZ.add(newMZs.get(index));
                                break;
                            default: //backbone and internal
                                StringBuilder sb = new StringBuilder();
                                ArrayList<Integer> letterPositions = new ArrayList<>();
                                ArrayList<Integer> numbers = new ArrayList<>();
                                letterPositions.add(ionName.length());
                                //get letters first
                                for (int i = ionName.length() - 2; i > -1; i--) { //last char will always be number
                                    char c = ionName.charAt(i);
                                    if (Character.isLetter(c)) {
                                        sb.append(c);
                                        letterPositions.add(i);
                                    }
                                }
                                //get numbers
                                for (int i = letterPositions.size() - 1; i > 0; i--) {
                                    numbers.add(Integer.parseInt(ionName.substring(letterPositions.get(i) + 1, letterPositions.get(i - 1))));
                                }
                                if (numbers.size() == 1) {
                                    newMZ.add(shiftedMC.calcMass(numbers.get(0), sb.toString(), charge, nl));
                                } else {
                                    newMZ.add(shiftedMC.calcMass(numbers.get(0), numbers.get(1), sb.toString(), nl));
                                }
                                break;
                        }
                    }

                    if (newMZ.size() == 1) {
                        finalMZs.add(newMZ.iterator().next());
                        finalIntensities.add(newIntensities.get(index));
                    }
                    index += 1;
                }

                //now can assign new PredictionEntry
                PredictionEntry newPred = new PredictionEntry();
                float[] mzArray = new float[finalMZs.size()];
                float[] intArray = new float[finalIntensities.size()];
                for (int i = 0; i < mzArray.length; i++) {
                    mzArray[i] = finalMZs.get(i);
                    intArray[i] = finalIntensities.get(i);
                }
                newPred.setMzs(mzArray);
                newPred.setIntensities(intArray);
                newPred.setRT(pe.RT);
                newPred.setIM(pe.IM);
                allPreds.put(lSplit[0] + "|" + lSplit[1], newPred);
            }
        }
    }
}
