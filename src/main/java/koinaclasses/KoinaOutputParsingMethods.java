package koinaclasses;

import java.util.List;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

//the annotation parsers of KoinaModelSpec.AnnotationParser. Each throws on an annotation it does not know, rather
//than adding a fragment with a missing or wrong ion type
public class KoinaOutputParsingMethods {
    private static final Pattern ION_NUMBER_CHARGE = Pattern.compile("([by])(\\d+)\\+(\\d+)");

    //"b3+1", "y12+2"
    public static void parseIonNumberCharge(String result, List<Integer> charges,
                                            List<String> fragmentIonTypes, List<Integer> fragNums) {
        Matcher m = ION_NUMBER_CHARGE.matcher(result);
        if (!m.matches()) {
            throw unknownAnnotation(result, "a b or y ion written as \"y3+1\"");
        }
        charges.add(Integer.parseInt(m.group(3)));
        fragmentIonTypes.add(m.group(1));
        fragNums.add(Integer.parseInt(m.group(2)));
    }

    public static void parseUnispec(String result, List<Integer> charges,
                                    List<String> fragmentIonTypes, List<Integer> fragNums) {
        //charge
        String[] info = result.split("\\^");
        if (info.length > 1) {
            charges.add(Integer.parseInt(info[1].substring(0, 1)));
        } else {
            charges.add(1);
        }

        //fragment ion type
        boolean NL = result.contains("-"); //neutral loss
        char first = result.charAt(0);
        switch (first) {
            case 'y':
                if (NL) {
                    fragmentIonTypes.add("y-NL");
                } else {
                    fragmentIonTypes.add("y");
                }
                break;
            case 'b':
                if (NL) {
                    fragmentIonTypes.add("b-NL");
                } else {
                    fragmentIonTypes.add("b");
                }
                break;
            case 'a':
                if (NL) {
                    fragmentIonTypes.add("a-NL");
                } else {
                    fragmentIonTypes.add("a");
                }
                break;

            case 'p':
                fragmentIonTypes.add("p");
                break;
            case 'I': //internal or immonium
                if (result.startsWith("Int")) {
                    if (NL) {
                        fragmentIonTypes.add("int-NL");
                    } else {
                        fragmentIonTypes.add("int");
                    }
                } else {
                    fragmentIonTypes.add("imm");
                }
                break;
            default:
                throw unknownAnnotation(result, "a UniSpec y, b, a, precursor, internal or immonium ion");
        }

        //fragment number
        switch (first) {
            case 'y':
            case 'b':
            case 'a':
                String num = "";
                for (int cidx = 1; cidx < result.length(); cidx++) {
                    char c = result.charAt(cidx);
                    if (Character.isDigit(c)) {
                        num += c;
                    } else {
                        break;
                    }
                }
                fragNums.add(Integer.parseInt(num));
                break;
            default:
                fragNums.add(0); //set as 0 if it's a more complicated fragment. Will need some full annotation here
                break;
        }
    }

    public static void parsePrositMultifrag(String result, List<Integer> charges,
                                            List<String> fragmentIonTypes, List<Integer> fragNums) {
        //charge
        String[] info = result.split("\\^");
        if (info.length > 1) {
            charges.add(Integer.parseInt(info[1].substring(0, 1)));
        } else {
            charges.add(1);
        }

        //fragment ion type
        char first = result.charAt(0);
        switch (first) {
            case 'A':
                fragmentIonTypes.add("a+1");
                break;
            case 'C':
                fragmentIonTypes.add("c");
                break;
            case 'X':
                fragmentIonTypes.add("x+1");
                break;
            case 'Z':
                fragmentIonTypes.add("zprime");
                break;
            case 'a':
                fragmentIonTypes.add("a");
                break;
            case 'b':
                fragmentIonTypes.add("b");
                break;
            case 'c':
                fragmentIonTypes.add("cdot");
                break;
            case 'x':
                fragmentIonTypes.add("x");
                break;
            case 'y':
                fragmentIonTypes.add("y");
                break;
            case 'z':
                fragmentIonTypes.add("zdot");
                break;
            default:
                throw unknownAnnotation(result, "a MultiFrag a, b, c, x, y or z ion");
        }

        //fragment number
        int number = Integer.parseInt(info[0].substring(1));
        fragNums.add(number);
    }

    //Altimeter's labels (mzPAF, Altimeter_2024_filter_splines/ion_dictionary.txt in Koina): a b or y ion and its
    //losses, then "^" and the charge if above 1 ("y5", "b3^2", "y5-H2O^2"); the precursor the same way ("p",
    //"p^3"); an immonium ion: "I", the residue, its mod in parentheses, and losses and gains ("IY",
    //"IC(Carbamidomethyl)", "IH-NH3+CO", "IH+2CO")
    private static final Pattern ALTIMETER_ION = Pattern.compile("([byp])(\\d*)((?:-[A-Za-z0-9]+)*)(?:\\^(\\d+))?");
    private static final Pattern ALTIMETER_IMMONIUM =
            Pattern.compile("I[A-Z](?:\\([A-Za-z]+\\))?(?:[+-]\\d*[A-Z][A-Za-z0-9]*)*");

    //losses are read as parseUnispec reads them: a b or y ion with one is b-NL or y-NL, the precursor is p and an
    //immonium ion imm with or without one
    public static void parseAltimeter(String result, List<Integer> charges,
                                      List<String> fragmentIonTypes, List<Integer> fragNums) {
        Matcher m = ALTIMETER_ION.matcher(result);
        if (m.matches()) {
            String type = m.group(1);
            boolean precursor = type.equals("p");
            //a b or y ion has a number, the precursor none
            if (precursor != m.group(2).isEmpty()) {
                throw unknownAnnotation(result, ALTIMETER_EXPECTED);
            }
            charges.add(m.group(4) == null ? 1 : Integer.parseInt(m.group(4)));
            if (precursor) {
                fragmentIonTypes.add("p");
                fragNums.add(0);
            } else {
                fragmentIonTypes.add(m.group(3).isEmpty() ? type : type + "-NL");
                fragNums.add(Integer.parseInt(m.group(2)));
            }
            return;
        }
        if (ALTIMETER_IMMONIUM.matcher(result).matches()) {
            charges.add(1);
            fragmentIonTypes.add("imm");
            fragNums.add(0);
            return;
        }
        throw unknownAnnotation(result, ALTIMETER_EXPECTED);
    }

    private static final String ALTIMETER_EXPECTED = "an Altimeter b, y, precursor or immonium ion written as "
            + "\"y3\", \"b5^2\", \"p^2\" or \"IY\"";

    private static IllegalArgumentException unknownAnnotation(String annotation, String expected) {
        return new IllegalArgumentException("Koina answered with the fragment annotation \"" + annotation
                + "\", which MSBooster cannot read; it expects " + expected);
    }
}
