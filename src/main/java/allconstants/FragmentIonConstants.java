package allconstants;

import features.spectra.MassCalculator;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

import static utils.Print.printError;

public class FragmentIonConstants {
    public static String ignoredFragmentIonTypes = ""; //split with commas
    public static String onlyFragmentIonTypes = ""; //split with commas
    public static String[] fragmentIonHierarchy;
    public static Set<String> fragmentIonHierarchySet;
    public static String divideFragments = "0";

    public static void makeFragmentIonHierarchy() {
        switch (Constants.FragmentationType) {
            case "HCD":
                fragmentIonHierarchy = new String[]{"p", "imm", "y", "b", "a", "p-NL",
                        "y-NL", "b-NL", "a-NL", "int", "int-NL", "unknown"};
            case "ETD":
                fragmentIonHierarchy = new String[]{"zdot", "c", "z", "y", "unknown"};
            case "ETHCD":
                fragmentIonHierarchy = new String[]{"imm", "y", "b", "a", "zdot", "c", "z", "cdot",
                        "y-NL", "b-NL", "a-NL", "int", "int-NL", "unknown"};
            default:  //everything else, like CID
                fragmentIonHierarchy = new String[]{"imm", "y", "b", "a",
                        "y-NL", "b-NL", "a-NL", "int", "int-NL", "unknown"};
        }
        makeLowestFragmentIonType();
        fragmentIonHierarchySet = new HashSet<>(Arrays.asList(fragmentIonHierarchy));
    }
    public static Set<String> makeIgnoredFragmentIonTypes() { //used here and also for experimental peaks in mgf format annotated with ion types
        Set<String> ignoredFragmentIonTypesSet = new HashSet<>();
        Set<String> onlyFragmentIonTypesSet = new HashSet<>();
        if (!onlyFragmentIonTypes.isEmpty()) {
            String[] commaSplit = onlyFragmentIonTypes.split(",");
            for (String s : commaSplit) {
                String fragmentIonType = s.trim();
                if (MassCalculator.allowedFragmentIonTypes.contains(fragmentIonType)) {
                    onlyFragmentIonTypesSet.add(fragmentIonType);
                } else {
                    printError(fragmentIonType + " is not a supported fragment ion type to include. " +
                            "Please choose from " + MassCalculator.allowedFragmentIonTypes);
                    System.exit(1);
                }
            }
            for (String fragment : MassCalculator.allowedFragmentIonTypes) {
                if (! onlyFragmentIonTypesSet.contains(fragment)) {
                    ignoredFragmentIonTypesSet.add(fragment);
                }
            }
        } else if (!ignoredFragmentIonTypes.isEmpty()) {
            //only filter if not excluding certain fragment ion types
            //check that this is allowed
            String[] commaSplit = ignoredFragmentIonTypes.split(",");
            for (String s : commaSplit) {
                String fragmentIonType = s.trim();
                if (MassCalculator.allowedFragmentIonTypes.contains(fragmentIonType)) {
                    ignoredFragmentIonTypesSet.add(fragmentIonType);
                } else {
                    printError(fragmentIonType + " is not a supported fragment ion type to exclude. " +
                            "Please choose from " + MassCalculator.allowedFragmentIonTypes);
                    System.exit(1);
                }
            }
        }
        return ignoredFragmentIonTypesSet;
    }
    private static void makeLowestFragmentIonType() { //use this to make sure fragment type is above priority threshold
        Set<String> ignoredFragmentIonTypesSet = makeIgnoredFragmentIonTypes();
        int index = 0;
        for (int i = fragmentIonHierarchy.length - 1; i > -1; i--) {
            String ion = fragmentIonHierarchy[i];
            if (! ignoredFragmentIonTypesSet.contains(ion)) {
                index = i;
                break;
            }
        }
        fragmentIonHierarchy = Arrays.copyOfRange(fragmentIonHierarchy, 0, index + 1);
    }
}
