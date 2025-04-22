package koinaclasses;

import java.util.ArrayList;

public class KoinaOutputParsingMethods {
    public static void parseUnispec(String result, ArrayList<String> fullAnnotation, ArrayList<Integer> charges,
                                    ArrayList<String> fragmentIonTypes, ArrayList<Integer> fragNums) {
        fullAnnotation.add(result);
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

    //TODO
    //Do not need full annotation
    public static void parsePrositMultifrag(String result, ArrayList<Integer> charges,
                                            ArrayList<String> fragmentIonTypes, ArrayList<Integer> fragNums) {
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
        }

        //fragment number
        int number = Integer.parseInt(info[0].substring(1));
        fragNums.add(number);
    }
}
