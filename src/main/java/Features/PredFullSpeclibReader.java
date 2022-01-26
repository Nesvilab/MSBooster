package Features;

import umich.ms.fileio.exceptions.FileParsingException;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;
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
//            //find locations of PTMs
//            ArrayList<Integer> starts = new ArrayList<>();
//            ArrayList<Integer> ends = new ArrayList<>();
//            for (int i = 0; i < l.length(); i++) {
//                if (l.charAt(i) == '[') {
//                    starts.add(i);
//                } else if (l.charAt(i) == ']') {
//                    ends.add(i);
//                }
//            }
//            boolean newEntryNeeded = false;
//            for (int i = 0; i < starts.size(); i++) {
//                Double reportedMass = Double.parseDouble(l.substring(starts.get(i) + 1, ends.get(i)));
//                if (!(reportedMass.equals(Constants.carbamidomethylationMass) || reportedMass.equals(Constants.oxidationMass))) {
//                    newEntryNeeded = true;
//                    break;
//                }
//            }
//
//            //given that C is assumed to be 57, need to check if no mod on it
//            if (! newEntryNeeded) {
//                for (int i = 0; i < l.length(); i++) {
//                    if (l.charAt(i) == 'C') {
//                        if (l.charAt(i + 1) != '[') {
//                            newEntryNeeded = true;
//                        }
//                    }
//                }
//            }
//
//            if (newEntryNeeded) {

            //doing fragment annotation for everything, not just modified ones
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

            //filter here so don't need to annotate everything
            double minIntensity = 0f;
            int numFragments = mzs.length;
            Set<String> ignoredFragmentIonTypes = new HashSet<>();

            if (Constants.ignoredFragmentIonTypes.equals("")) {
                if (Constants.useBasePeak) {
                    for (float i : intensities) {
                        if (i > minIntensity) {
                            minIntensity = i;
                        }
                    }
                    minIntensity = minIntensity * Constants.percentBasePeak / 100f;
                }
                if (Constants.useTopFragments) {
                    numFragments = Math.min(numFragments, Constants.topFragments);
                }
            } else { //only filter if not excluding certain fragment ion types
                //check that this is allowed
                String[] commaSplit = Constants.ignoredFragmentIonTypes.split(",");
                for (int i = 0; i < commaSplit.length; i++) {
                    String fragmentIonType = commaSplit[i].trim();
                    if (MassCalculator.allowedFragmentIonTypes.contains(fragmentIonType)) {
                        ignoredFragmentIonTypes.add(fragmentIonType);
                    } else {
                        System.out.println(fragmentIonType + " is not a supported fragment ion type to exclude. " +
                                "Please choose from " + MassCalculator.allowedFragmentIonTypes);
                        System.exit(-1);
                    }
                }
            }

            //for each predicted peak, see if it can be annotated
            ArrayList<Float> newMZs = new ArrayList<>();
            ArrayList<Float> newIntensities = new ArrayList<>();
            ArrayList<Float> finalMZs = new ArrayList<>();
            ArrayList<Float> finalIntensities = new ArrayList<>();

            for (int i = 0; i < mzs.length; i++) {
                if (intensities[i] >= minIntensity) {
                    newIntensities.add(intensities[i]);
                    newMZs.add(mzs[i]);
                }
            }
            if (numFragments < newMZs.size()) {
                ArrayList<Float> newnewMZs = new ArrayList<>();
                ArrayList<Float> newnewIntensities = new ArrayList<>();
                for (int i = 0; i < numFragments; i++) {
                    int index = newIntensities.indexOf(Collections.max(newIntensities));
                    newnewMZs.add(newMZs.get(index));
                    newnewIntensities.add(newIntensities.get(index));
                    newIntensities.set(index, -1f);
                }
                newMZs = newnewMZs;
                newIntensities = newnewIntensities;
            }

            //annotate ion
            //any calculated annotated fragment is considered
            String[][] info = mc.annotateMZs(newMZs);
            String[] annotations = info[0];
            String[] fragmentIonTypes = info[1];
            //System.out.println(Arrays.toString(annotations));

            //for annotated fragment ions, shift their predicted peak
            //shift only if including the amino acid
            //create new massCalculator object, and change m/z to new one
            //give String in annotation to calcMass and round to closest tenth of m/z
            //adjust in newMZs
            MassCalculator shiftedMC = new MassCalculator(lSplit[0], lSplit[1]);
            int index = 0;
            for (int j = 0; j < annotations.length; j++) {
                //skip fragment ion if of ignored type
                if (ignoredFragmentIonTypes.contains(fragmentIonTypes[j])) {
                    continue;
                }

                String fragment = annotations[j];
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
                    String[] nlSplitSplit = nlSplit[0].split(":");
                    String ionName = nlSplitSplit[0];
                    switch (ionName) {
                        case "MH":
                            newMZ.add(shiftedMC.calcMass(shiftedMC.modMasses.size() - 2, "y", charge, nl));
                            break;
                        case "immonium": //TODO: why is immonium showing up twice?
                            //need to check if mox of modified and unmodified, and exclude if so
                            //don't know if both mod and unmod m/z will be that intensity, safer to exclude
                            Set<Double> mods = new HashSet<>();
                            char aa = nlSplitSplit[1].charAt(0);
                            for (int i = 0; i < shiftedMC.peptide.length(); i++) {
                                if (shiftedMC.peptide.charAt(i) == aa) {
                                    mods.add(shiftedMC.modMasses.get(i + 1));
                                }
                            }
                            for (int i = 0; i < mods.size(); i++) {
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
                if (newMZ.size() == 1) { //problem if one possibility is shifted and other isn't
                    finalMZs.add(newMZ.iterator().next());
                    finalIntensities.add(newIntensities.get(index));
                } else {
                    float totalFloat = 0f;
                    float compareFloat1 = newMZ.iterator().next();
                    totalFloat += compareFloat1;
                    boolean problem = false;
                    for (int i = 0; i < newMZ.size() - 1; i++) {
                        float compareFloat2 = newMZ.iterator().next();
                        totalFloat += compareFloat2;
                        if (Math.abs(compareFloat1 - compareFloat2) > 0.1f) {
                            problem = true;
                        }
                    }
                    if (! problem) {
                        finalMZs.add(totalFloat / newMZ.size());
                        finalIntensities.add(newIntensities.get(index));
                    }
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
            newPred.setFragmentIonTypes(fragmentIonTypes);
            newPred.setMassCalculator(shiftedMC);
            allPreds.put(lSplit[0] + "|" + lSplit[1], newPred);
//            }
        }
    }
}
