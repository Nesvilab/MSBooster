package predictions;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.ArrayList;

import static features.spectra.MassCalculator.allNeutralLossMasses;
import static utils.NumericUtils.isUppercaseOrDigit;

//class to parse annotations for unispec and predfull
//designed to handle output from unispec koina prediction and read in dictionary of unispec allowed fragments
//TODO: support for c/z ions
public class FragmentAnnotationParser {
    public String fullAnnotation;
    public int charge = 1;
    public String fragmentIonType;
    public int fragnum = 0;
    public int internalStartPosition;
    public int internalExtent;
    public int isotope = 0;
    public String internalSequence;
    public float neutralLoss = 0f;
    public ArrayList<String> neutralLossStrings = new ArrayList<>();

    public FragmentAnnotationParser(String annotation) throws IOException, URISyntaxException {
        annotation = annotation.split(";")[0]; //if multiple possible annotations (from default fragment assignment), use first
        fullAnnotation = annotation;

        //charge
        String[] info = annotation.split("\\^");
        if (info.length > 1) {
            charge = Integer.parseInt(info[1].substring(0, 1));
        }

        //fragment ion type
        boolean NL = annotation.contains("-"); //neutral loss
        char first = annotation.charAt(0);
        switch (first) {
            case 'y':
                if (NL) {
                    fragmentIonType = "y-NL";
                } else {
                    fragmentIonType = "y";
                }
                break;
            case 'b':
                if (NL) {
                    fragmentIonType = "b-NL";
                } else {
                    fragmentIonType = "b";
                }
                break;
            case 'a':
                if (NL) {
                    fragmentIonType = "a-NL";
                } else {
                    fragmentIonType = "a";
                }
                break;
            case 'p':
                if (! annotation.contains("^")) {
                    charge = 0;
                }
                if (NL) {
                    fragmentIonType = "p-NL";
                } else {
                    fragmentIonType = "p";
                }
                break;
            case 'I': //internal or immonium
                if (annotation.startsWith("Int")) {
                    if (NL) {
                        fragmentIonType = "int-NL";
                    } else {
                        fragmentIonType = "int";
                    }
                } else {
                    fragmentIonType = "imm";
                }
                break;
            case 'u': //unknown
                fragmentIonType = "unknown";
                break;
            default:
                utils.Print.printError(annotation + " not supported for parsing. Exiting");
                System.exit(1);
        }

        //fragment number
        switch (first) {
            case 'y':
            case 'b':
            case 'a':
                String num = "";
                for (int cidx = 1; cidx < annotation.length(); cidx++) {
                    char c = annotation.charAt(cidx);
                    if (Character.isDigit(c)) {
                        num += c;
                    } else {
                        break;
                    }
                }
                fragnum = Integer.parseInt(num);
                break;
        }
        if (annotation.startsWith("Int")) {
            //start>extent
            //in MassCalculator
            //num1 is peptide length - start
            //num2 should be extent
            if (annotation.contains(":")) {
                //formatting from default possible fragments method
            } else if (annotation.contains("/")) { //koina output
                String[] annotationSplit = annotation.split("/");
                internalSequence = annotationSplit[1];
                internalStartPosition = Integer.parseInt(annotationSplit[2]);

                int plusIndex = internalSequence.indexOf('+');
                int minusIndex = internalSequence.indexOf('-');
                int firstIndex = (plusIndex == -1) ? minusIndex : (minusIndex == -1) ? plusIndex : Math.min(plusIndex, minusIndex);
                if (firstIndex == -1) {
                    internalExtent = internalSequence.length();
                } else {
                    internalExtent = firstIndex;
                }
            } else { //text based fragment dictionary file
                String internalString = annotation.substring(3);
                String[] internalStringSplit = internalString.split(">");
                internalStartPosition = Integer.parseInt(internalStringSplit[0]);
                String internalExtentString = "";
                for (int i = 0; i < internalStringSplit[1].length(); i++) {
                    char currentChar = internalStringSplit[1].charAt(i);
                    if (Character.isDigit(currentChar)) {
                        internalExtentString += currentChar;
                    } else {
                        break;
                    }
                }
                internalExtent = Integer.parseInt(internalExtentString);
            }
        }

        //isotope
        String isotopeString = annotation;
        if (internalSequence != null) {
            isotopeString = internalSequence;
        }
        if (isotopeString.endsWith("i")) {
            char num = isotopeString.charAt(isotopeString.length() - 2);
            if (num == '+') {
                isotope = 1;
            } else {
                isotope = Integer.parseInt(String.valueOf(num));
            }
        }

        //neutral loss
        //allowed characters ae uppercase letter, numbers, +-
        String breakChars = "i^/";
        String[] annotationSplit = annotation.split("-");
        if (annotationSplit.length != 1) {
            ArrayList<Integer> multipliers = new ArrayList<>();

            for (int split = 1; split < annotationSplit.length; split++) {
                String nlString = annotationSplit[split];
                String nl = "";
                int multiplier = 1;

                if (Character.isDigit(nlString.charAt(0))) {
                    multiplier *= Integer.parseInt(String.valueOf(nlString.charAt(0)));
                } else {
                    nl += nlString.charAt(0);
                }
                for (int charI = 1; charI < nlString.length(); charI++) {
                    char currentChar = nlString.charAt(charI);
                    if (isUppercaseOrDigit(currentChar)) {
                        nl += currentChar;
                    } else if (breakChars.indexOf(currentChar) != -1) {
                        break;
                    } else if (currentChar == '+') {
                        if (! nl.isEmpty()) {
                            neutralLossStrings.add(nl);
                            multipliers.add(multiplier);
                        }
                        nl = "";
                        multiplier = -1;
                        if (Character.isDigit(nlString.charAt(charI + 1))) {
                            multiplier *= Integer.parseInt(String.valueOf(nlString.charAt(charI + 1)));
                            charI++;
                        }
                    } else {
                        utils.Print.printError(currentChar + " not recognized. Exiting");
                        System.exit(1);
                    }
                }
                if (! nl.isEmpty()) {
                    neutralLossStrings.add(nl);
                    multipliers.add(multiplier);
                }
            }

            //query neutral losses and calculate mass
            for (int i = 0; i < neutralLossStrings.size(); i++) {
                neutralLoss += allNeutralLossMasses.get(neutralLossStrings.get(i)) * multipliers.get(i);
            }
        }
    }

    public static void main(String[] args) {
        String fileName = "Z:/yangkl/manuscripts/ImmunopeptidomeMethods/UnispecVsPredfull/unispec_fragments.txt";
        //String fileName = "Z:/yangkl/manuscripts/ImmunopeptidomeMethods/UnispecVsPredfull/koina_output_example.txt";

        try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {
            String line;
            while ((line = br.readLine()) != null) {
                line = line.replace("\"", "");
                line = line.replace(",", "");
                new FragmentAnnotationParser(line);
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (URISyntaxException e) {
            throw new RuntimeException(e);
        }
    }
}
