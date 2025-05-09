package allconstants;

import predictions.FragmentAnnotationParser;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Objects;
import java.util.Set;
import java.util.TreeSet;

import static utils.Print.printError;
import static utils.Print.printInfo;

public class FragmentIonConstants implements ConstantsInterface {
    public static final HashSet<String> allowedFragmentIonTypes = new HashSet<>(Arrays.asList(
            "z", "c", "y", "a", "x", "b", "zdot", "cdot", "zprime", "x+1", "a+1",
            "p", "imm", "int", "int-NL",
            "z-NL", "c-NL", "y-NL", "a-NL", "x-NL", "b-NL",
            "p-NL", "unknown"));

    //these parameters handle the fragments used, with decreasing priority
    public static String ignoredFragmentIonTypes = ""; //split with commas
    public static String onlyFragmentIonTypes = ""; //split with commas
    public static String[] fragmentIonHierarchy; //ordering fragments in decreasing priority
    public static HashSet<String> fragmentIonHierarchySet;

    //this is how we group fragments to similarity calculations
    public static Integer divideFragments = 0;
    public static String customFragmentGroups;
    public static TreeSet<String>[] fragmentGroups = new TreeSet[1]; //iterate through sets
    public static TreeSet<String>[] createFragmentGroups() {
        TreeSet[] fg;
        switch (divideFragments) {
            case 0: //use everything in one score
                fg = new TreeSet[1];
                fg[0] = new TreeSet<>();
                for (String s : fragmentIonHierarchySet) {
                    fg[0].add(s);
                }
                break;
            case 1: //primary & aux separation
                fg = new TreeSet[2];
                fg[0] = new TreeSet<>();
                fg[1] = new TreeSet<>();
                for (String s : primaryFragmentIonTypes) {
                    fg[0].add(s);
                }
                for (String s : auxFragmentIonTypes) {
                    fg[1].add(s);
                }
                break;
            case 2: //break into y/b and others, HCD/CID version
                fg = new TreeSet[2];
                fg[0] = new TreeSet<>();
                fg[1] = new TreeSet<>();
                fg[0].add("y");
                fg[0].add("b");
                for (String s : fragmentIonHierarchy) {
                    if (!Objects.equals(s, "y") && !Objects.equals(s, "b")) {
                        fg[1].add(s);
                    }
                }
                break;
            case 3: //break into y/b, others, and unknown
                fg = new TreeSet[3];
                fg[0] = new TreeSet<>();
                fg[1] = new TreeSet<>();
                fg[2] = new TreeSet<>();
                fg[0].add("y");
                fg[0].add("b");
                for (String s : fragmentIonHierarchy) {
                    if (!Objects.equals(s, "y") && !Objects.equals(s, "b") && !Objects.equals(s, "unknown")) {
                        fg[1].add(s);
                    }
                }
                fg[2].add("unknown");
                break;
            case 4: //all individual groups
                fg = new TreeSet[fragmentIonHierarchy.length];
                int j = 0;
                for (String s : fragmentIonHierarchy) {
                    fg[j] = new TreeSet<>();
                    fg[j].add(s);
                    j++;
                }
                break;
            case 5: //break into zdot/z/c/y and unknown, ETD version
                fg = new TreeSet[2];
                fg[0] = new TreeSet<>();
                fg[1] = new TreeSet<>();
                fg[0].add("zdot");
                fg[0].add("z");
                fg[0].add("c");
                fg[0].add("y");
                fg[1].add("unknown");
                break;
            default: //custom
                printError("Custom fragment groups not yet supported. Exiting");
                System.exit(1);
                String[] splits = customFragmentGroups.split(";");
                fg = new TreeSet[splits.length];
                for (int i = 0; i < splits.length; i++) {
                    fg[i] = new TreeSet<>();
                    for (String s : splits[i].split(",")) {
                        fg[i].add(s);
                    }
                }
                break;
        }
        return fg;
    }

//            if (Constants.adaptiveFragmentNum) {
//                Constants.topFragments = 36; //TODO think of better way than hardcoding
//            } else if (Constants.divideFragments.equals("1")) { //standard setting of yb vs others
//                Constants.divideFragments = "y_b;immonium_a_y-NL_b-NL_a-NL_internal_internal-NL_unknown";
//                Constants.topFragments = 12;
//            } else if (Constants.divideFragments.equals("2")) {
//                Constants.divideFragments = "y;b;immonium;a;y-NL;b-NL;a-NL;internal;internal-NL;unknown";
//                Constants.topFragments = 6;
//            } else if (Constants.divideFragments.equals("3")) { //standard setting of yb vs others
//                Constants.divideFragments = "y_b_y-NL_b-NL;immonium_a_a-NL_internal_internal-NL_unknown";
//                Constants.topFragments = 12;
//            } else if (Constants.divideFragments.equals("4")) { //etd
//                Constants.divideFragments = "c_z;zdot_y_unknown";
//                Constants.topFragments = 12;
//            } else if (Constants.divideFragments.equals("5")) { //ethcd
//                Constants.divideFragments = "b_y_c_z;immonium_a_cdot_zdot_y-NL_b-NL_a-NL_internal_internal-NL_unknown";
//                Constants.topFragments = 12;
//            } else if (Constants.divideFragments.equals("0") && Constants.spectraModel.equals("DIA-NN")) {
//                Constants.topFragments = 20;
//            }

    //TODO: allow defining fragment ion hierarchy
    public static final HashSet<String> allowedFragmentationTypes = new HashSet<>(Arrays.asList(
            "HCD", "CID", "ETD", "ETHCD", "ECD", "EID", "UVPD", "ETCID"));
    public static void makeFragmentIonHierarchy() {
        switch (Constants.FragmentationType.toUpperCase()) {
            case "HCD":
            case "CID":
                fragmentIonHierarchy = new String[]{"imm", "y", "b", "a",
                        "y-NL", "b-NL", "a-NL", "int", "int-NL", "unknown"};
                break;
            case "ETD":
                fragmentIonHierarchy = new String[]{"zdot", "c", "z", "y", "unknown"};
                break;
            case "ETHCD":
                fragmentIonHierarchy = new String[]{"imm", "y", "b", "a", "zdot", "c", "z", "cdot",
                        "y-NL", "b-NL", "a-NL", "int", "int-NL", "unknown"};
                break;
            case "ECD":
                fragmentIonHierarchy = new String[]{"c", "zdot", "zprime", "y", "cdot", "b", "a+1",
                        "c", "x+1", "a", "x", "unknown"};
                break;
            case "EID":
                fragmentIonHierarchy = new String[]{"y", "b", "a", "zdot", "c", "a+1", "x+1",
                        "zprime", "cdot", "x", "unknown"};
                break;
            case "UVPD":
                fragmentIonHierarchy = new String[]{"y", "b", "a", "a+1", "zdot", "c",
                        "x+1", "zprime", "cdot", "x", "unknown"};
                break;
            case "ETCID":
                fragmentIonHierarchy = new String[]{"cdot", "zdot", "zprime", "y", "b",
                        "a+1", "c", "x+1", "a", "x", "unknown"};
                break;
            default:  //everything else
                printInfo(Constants.FragmentationType + " not recognized. Setting fragment ion hierarchy to default.");
                fragmentIonHierarchy = new String[]{"imm", "y", "b", "a",
                        "y-NL", "b-NL", "a-NL", "int", "int-NL", "unknown"};
                break;
        }
        removeIgnoredFragments();
        fragmentIonHierarchySet = new HashSet<>(Arrays.asList(fragmentIonHierarchy));
        auxFragmentIonTypes.addAll(fragmentIonHierarchySet);
    }
    public static Set<String> makeIgnoredFragmentIonTypes() {
        Set<String> ignoredFragmentIonTypesSet = new HashSet<>();
        Set<String> onlyFragmentIonTypesSet = new HashSet<>();
        if (!onlyFragmentIonTypes.isEmpty()) {
            String[] commaSplit = onlyFragmentIonTypes.split(",");
            for (String s : commaSplit) {
                String fragmentIonType = s.trim();
                if (allowedFragmentIonTypes.contains(fragmentIonType)) {
                    onlyFragmentIonTypesSet.add(fragmentIonType);
                } else {
                    printError(fragmentIonType + " is not a supported fragment ion type to include. " +
                            "Please choose from " + allowedFragmentIonTypes);
                    System.exit(1);
                }
            }
            for (String fragment : allowedFragmentIonTypes) {
                if (! onlyFragmentIonTypesSet.contains(fragment)) {
                    ignoredFragmentIonTypesSet.add(fragment);
                }
            }
        }
        if (!ignoredFragmentIonTypes.isEmpty()) {
            //only filter if not excluding certain fragment ion types
            //check that this is allowed
            String[] commaSplit = ignoredFragmentIonTypes.split(",");
            for (String s : commaSplit) {
                String fragmentIonType = s.trim();
                if (allowedFragmentIonTypes.contains(fragmentIonType)) {
                    ignoredFragmentIonTypesSet.add(fragmentIonType);
                } else {
                    printError(fragmentIonType + " is not a supported fragment ion type to exclude. " +
                            "Please choose from " + allowedFragmentIonTypes);
                    System.exit(1);
                }
            }
        }
        return ignoredFragmentIonTypesSet;
    }
    private static void removeIgnoredFragments() { //use this to make sure fragment type is above priority threshold
        Set<String> ignoredFragmentIonTypesSet = makeIgnoredFragmentIonTypes();
        ArrayList<String> fragments = new ArrayList<>();
        for (String ion : fragmentIonHierarchy) {
            if (!ignoredFragmentIonTypesSet.contains(ion)) {
                fragments.add(ion);
            }
        }
        fragmentIonHierarchy = new String[fragments.size()];
        for (int i = 0; i < fragments.size(); i++) {
            fragmentIonHierarchy[i] = fragments.get(i);
        }
    }

    //contains fragment ion types predicted by primary, not auxiliary spectra model
    public static HashSet<String> primaryFragmentIonTypes = new HashSet<>();
    public static HashSet<String> auxFragmentIonTypes = new HashSet<>();

    public static void setPrimaryAndAuxFragmentIonTypes(String[] primaryTypes) {
        List<String> primaryTypesList = Arrays.asList(primaryTypes);
        primaryFragmentIonTypes.addAll(primaryTypesList);
        primaryTypesList.forEach(auxFragmentIonTypes::remove);
    }

    //to be used by possible unispec mzs method
    public static ArrayList<FragmentAnnotationParser> fragmentAnnotationParserArrayList;

    static {
        try {
            fragmentAnnotationParserArrayList = makeFragmentAnnotationParserArrayList();
        } catch (IOException | URISyntaxException e) {
            throw new RuntimeException(e);
        }
    }

    public static ArrayList<FragmentAnnotationParser> makeFragmentAnnotationParserArrayList() throws IOException, URISyntaxException {
        final InputStream stream = FragmentIonConstants.class.getClassLoader().getResourceAsStream(
                "fragment_annotation/unispec_fragments.txt");
        final InputStreamReader reader = new InputStreamReader(stream);
        final BufferedReader fragmentsFile = new BufferedReader(reader);
        ArrayList<FragmentAnnotationParser> faps = new ArrayList<>();

        String line;
        while((line = fragmentsFile.readLine()) != null) {
            FragmentAnnotationParser fap = new FragmentAnnotationParser(line);
            faps.add(fap);
        }
        fragmentsFile.close();

        return faps;
    }

    public static Boolean annotatePredfullLikeUnispec = true;
}
