/*
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
 * implied. See the License for the specific language governing
 * permissions and limitations under the License.
 */

package features.spectra;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import allconstants.FragmentIonConstants;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeSet;
import mainsteps.MzmlScanNumber;
import mainsteps.PeptideObj;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import predictions.FragmentAnnotationParser;
import utils.StatMethods;

class SpectrumComparisonTest {

    private static final float[] EXP_MZS = {100f, 200f, 300f, 400f};
    private static final float[] PRED_INTENSITIES = {1f, 0.8f, 0.6f, 0.4f};

    @BeforeEach
    void setUp() {
        //PeptideObj's static baseMap iterates this; populate before PeptideObj loads
        FragmentIonConstants.fragmentIonHierarchy = new String[]{"y", "b"};
        FragmentIonConstants.primaryFragmentIonTypes.clear(); //-> matchWithDaltonsDefault decides
        Constants.matchWithDaltonsDefault = false;            //match on ppm
        Constants.ppmTolerance = 20f / 1000000f;
        Constants.minFragments = 2;
    }

    private SpectrumComparison compare(float[] expIntensities, float[] predMZs) {
        MzmlScanNumber scan = new MzmlScanNumber(1, EXP_MZS.clone(), expIntensities.clone(), 1.0f, 0f);
        PeptideObj pepObj = new PeptideObj();
        pepObj.scanNumObj = scan;
        return new SpectrumComparison(pepObj, scan.getExpMZs(), scan.getExpIntensities(),
                predMZs, PRED_INTENSITIES.clone(), new String[]{"y", "y", "y", "y"}, false);
    }

    // top6matchedIntensity() normalizes by the summed log intensity of every peak STRICTLY above
    // the spectrum mean. When all peaks share one intensity (a SCIEX TOF quantization artifact)
    // no peak clears the mean, so that denominator is 0 and the division escaped as Infinity/NaN
    // into the pin, breaking Percolator and crashing the xchart score histograms.
    @Test
    void flatSpectrumWithMatchesIsNotInfinity() {
        //every predicted fragment matches -> score > 0 -> score/0 was Infinity
        SpectrumComparison sc = compare(new float[]{48f, 48f, 48f, 48f}, new float[]{100f, 200f, 300f, 400f});

        double score = sc.top6matchedIntensity();

        assertTrue(Double.isFinite(score), "flat spectrum must not yield a non-finite score, was " + score);
        assertEquals(0.0, score, "an unnormalizable spectrum should fall back to 0, like the minFragments guard");
    }

    @Test
    void flatSpectrumWithoutMatchesIsNotNaN() {
        //no predicted fragment matches -> score == 0 -> 0.0/0.0 was NaN
        SpectrumComparison sc = compare(new float[]{48f, 48f, 48f, 48f}, new float[]{111f, 222f, 333f, 444f});

        double score = sc.top6matchedIntensity();

        assertTrue(Double.isFinite(score), "flat spectrum must not yield a non-finite score, was " + score);
        assertEquals(0.0, score, "an unnormalizable spectrum should fall back to 0, like the minFragments guard");
    }

    @Test
    void emptySpectrumIsNotNaN() {
        //mean() of an empty array is 0/0 = NaN, so no peak clears the threshold either
        MzmlScanNumber scan = new MzmlScanNumber(1, new float[0], new float[0], 1.0f, 0f);
        PeptideObj pepObj = new PeptideObj();
        pepObj.scanNumObj = scan;
        SpectrumComparison sc = new SpectrumComparison(pepObj, new float[0], new float[0],
                EXP_MZS.clone(), PRED_INTENSITIES.clone(), new String[]{"y", "y", "y", "y"}, false);

        double score = sc.top6matchedIntensity();

        assertTrue(Double.isFinite(score), "empty spectrum must not yield a non-finite score, was " + score);
        assertEquals(0.0, score);
    }

    // guard against the fix flattening ordinary spectra: one peak clears the mean here, so the
    // feature must still normalize as before rather than short-circuit to 0
    @Test
    void ordinarySpectrumStillNormalizes() {
        SpectrumComparison sc = compare(new float[]{10f, 20f, 30f, 100f}, new float[]{100f, 200f, 300f, 400f});

        double score = sc.top6matchedIntensity();

        //denominator = log(101); score = log(11)+log(21)+log(31)+log(101)
        double expected = (Math.log(11) + Math.log(21) + Math.log(31) + Math.log(101))
                / (float) Math.log(101);
        assertEquals(expected, score, 1e-6, "ordinary spectra must keep their existing normalization");
    }

    // bootstrapSimilarity draws the fragments of every bootstrap from a Random seeded with the scan
    // number, so that its score does not depend on the threads: one seed must always draw the same ones
    @Test
    void pickedPredictedDrawsTheSameFragmentsForOneSeed() {
        SpectrumComparison sc = compare(new float[]{10f, 20f, 30f, 100f}, new float[]{100f, 200f, 300f, 400f});

        for (long seed = 0; seed < 20; seed++) {
            SpectrumComparison first = sc.pickedPredicted(new Random(seed));
            SpectrumComparison second = sc.pickedPredicted(new Random(seed));

            assertTrue(first.predIntensities.length > 0, "seed " + seed + " drew no fragments");
            assertArrayEquals(first.predIntensities, second.predIntensities, "seed " + seed);
            assertArrayEquals(first.matchedIntensities, second.matchedIntensities, "seed " + seed);
        }
    }

    // bootstrapSimilarity must be one median over all Constants.bootstraps draws from one Random seeded with
    // the scan number, whatever Constants.numThreads is. It used to be the median of one thread's share
    @Test
    void bootstrapSimilarityIsOneMedianOverAllBootstraps() {
        float[] expMZs = {100f, 200f, 300f, 400f, 500f, 600f, 700f, 800f};
        float[] expIntensities = {10f, 50f, 20f, 80f, 5f, 40f, 30f, 60f};
        MzmlScanNumber scan = new MzmlScanNumber(1, expMZs, expIntensities, 1.0f, 0f);
        PeptideObj pepObj = new PeptideObj();
        pepObj.scanNumObj = scan;
        String[] fragmentIonTypes = new String[10];
        Arrays.fill(fragmentIonTypes, "y");
        SpectrumComparison sc = new SpectrumComparison(pepObj, scan.getExpMZs(), scan.getExpIntensities(),
                new float[]{100f, 200f, 300f, 400f, 500f, 600f, 700f, 800f, 900f, 1000f},
                new float[]{1f, 0.3f, 0.7f, 0.2f, 0.9f, 0.5f, 0.1f, 0.6f, 0.8f, 0.4f}, fragmentIonTypes, false);
        long seed = 1234;

        Random random = new Random(seed);
        ArrayList<Double> scores = new ArrayList<>();
        for (int i = 0; i < Constants.bootstraps; i++) {
            scores.add(sc.pickedPredicted(random).unweightedSpectralEntropy());
        }
        double expected = StatMethods.medianDouble(new ArrayList<>(scores));
        double oneShareOfEight = StatMethods.medianDouble(new ArrayList<>(scores.subList(0, Constants.bootstraps / 8)));
        assertNotEquals(expected, oneShareOfEight, "the spectrum must tell one thread's share of 8 from all bootstraps");

        Integer numThreads = Constants.numThreads;
        try {
            for (int threads : new int[]{1, 8}) {
                Constants.numThreads = threads;
                assertEquals(expected, sc.bootstrapSimilarity(seed), "numThreads " + threads);
            }
        } finally {
            Constants.numThreads = numThreads;
        }
    }

    // the fragments a spectrum could have (hypergeometricProbability and intersection) come from UniSpec's fragment
    // list for the ion types a model that annotates its fragments as UniSpec does predicted, as the primary or the
    // auxiliary model. For the types any other model predicted they are calculated, Altimeter's included, though its
    // m/z are Koina's too. annotatePredfullLikeUnispec takes UniSpec's list for every model
    @Test
    void thePossibleFragmentsAreUniSpecsForAModelThatAnnotatesAsUniSpecDoes() throws Exception {
        boolean annotateLikeUnispec = FragmentIonConstants.annotatePredfullLikeUnispec;
        ArrayList<FragmentAnnotationParser> unispecFragments = FragmentIonConstants.fragmentAnnotationParserArrayList;
        Set<String> primaryTypes = new HashSet<>(FragmentIonConstants.primaryFragmentIonTypes);
        HashSet<String> hierarchy = FragmentIonConstants.fragmentIonHierarchySet;
        Boolean matchWithDaltons = Constants.matchWithDaltons;
        Boolean matchWithDaltonsAux = Constants.matchWithDaltonsAux;
        String spectraModel = Constants.spectraModel;
        String auxSpectraModel = Constants.auxSpectraModel;
        String fragmentationType = Constants.FragmentationType;
        try {
            FragmentIonConstants.annotatePredfullLikeUnispec = false;
            FragmentIonConstants.fragmentAnnotationParserArrayList = new ArrayList<>(List.of(
                    new FragmentAnnotationParser("y1"), new FragmentAnnotationParser("y2"),
                    new FragmentAnnotationParser("b2")));
            FragmentIonConstants.primaryFragmentIonTypes.add("y"); //b ions come from the auxiliary model
            FragmentIonConstants.fragmentIonHierarchySet = new HashSet<>(Set.of("y", "b"));
            Constants.matchWithDaltons = false;
            Constants.matchWithDaltonsAux = false;
            Constants.FragmentationType = "HCD";

            MassCalculator mc = new MassCalculator("PEPTIDEK", 2);
            int calculatedY = new TreeSet<>(mc.possibleFragmentIons("y")).size();
            TreeSet<Float> calculatedYB = new TreeSet<>(mc.possibleFragmentIons("y"));
            calculatedYB.addAll(mc.possibleFragmentIons("b"));
            TreeSet<Float> calculatedYUnispecB = new TreeSet<>(mc.possibleFragmentIons("y"));
            calculatedYUnispecB.addAll(mc.possibleUnispecMzs("b"));
            assertTrue(calculatedY > 2, "y1 to y7 at charges 1 and 2, and their isotopes: " + calculatedY);

            assertEquals(2, possibleFragments("UniSpec", "", "y"), "y1 and y2");
            assertEquals(calculatedY, possibleFragments("Altimeter_2024_intensities", "", "y"));
            assertEquals(calculatedY, possibleFragments("Prosit_2020_intensity_HCD", "", "y"));
            assertEquals(calculatedYUnispecB.size(),
                    possibleFragments("Prosit_2020_intensity_HCD", "UniSpec", "y", "b"));
            assertEquals(calculatedYB.size(), possibleFragments("Prosit_2020_intensity_HCD",
                    "Altimeter_2024_intensities", "y", "b"));
            assertEquals(3, possibleFragments("UniSpec", "UniSpec", "y", "b"));

            FragmentIonConstants.annotatePredfullLikeUnispec = true;
            assertEquals(2, possibleFragments("Altimeter_2024_intensities", "", "y"));
        } finally {
            FragmentIonConstants.annotatePredfullLikeUnispec = annotateLikeUnispec;
            FragmentIonConstants.fragmentAnnotationParserArrayList = unispecFragments;
            FragmentIonConstants.primaryFragmentIonTypes.clear();
            FragmentIonConstants.primaryFragmentIonTypes.addAll(primaryTypes);
            FragmentIonConstants.fragmentIonHierarchySet = hierarchy;
            Constants.matchWithDaltons = matchWithDaltons;
            Constants.matchWithDaltonsAux = matchWithDaltonsAux;
            Constants.spectraModel = spectraModel;
            Constants.auxSpectraModel = auxSpectraModel;
            Constants.FragmentationType = fragmentationType;
        }
    }

    //how many fragments PEPTIDEK at charge 2 could have, for a prediction of fragments of these types by these models
    private static int possibleFragments(String spectraModel, String auxSpectraModel, String... types) {
        Constants.spectraModel = spectraModel;
        Constants.auxSpectraModel = auxSpectraModel;
        MzmlScanNumber scan = new MzmlScanNumber(1, EXP_MZS.clone(), new float[]{10f, 20f, 30f, 100f}, 1.0f, 0f);
        scan.ms2LowerLimit = 0.0;
        scan.ms2UpperLimit = 5000.0;
        PeptideObj pepObj = new PeptideObj();
        pepObj.scanNumObj = scan;
        pepObj.name = "PEPTIDEK|2";
        pepObj.charge = 2;
        String[] predTypes = new String[EXP_MZS.length];
        for (int i = 0; i < predTypes.length; i++) {
            predTypes[i] = types[i % types.length];
        }
        SpectrumComparison sc = new SpectrumComparison(pepObj, scan.getExpMZs(), scan.getExpIntensities(),
                new float[]{111f, 222f, 333f, 444f}, PRED_INTENSITIES.clone(), predTypes, false);
        sc.hypergeometricProbability();
        return sc.allMatchedIntensities.length;
    }
}
