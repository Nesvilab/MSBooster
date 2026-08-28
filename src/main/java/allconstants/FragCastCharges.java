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

package allconstants;

import java.util.concurrent.atomic.AtomicInteger;

import static utils.Print.printInfo;

/**
 * Which precursor charges FragCast can be asked to predict.
 *
 * <p>FragCast's peptide list takes an optional charge per row: a row that names one becomes exactly
 * that precursor, and a row that leaves it blank is expanded across {@code --min-charge}..
 * {@code --max-charge}. MSBooster always knows the charge - it comes from the pin - so it always
 * names it, and the predicted library holds exactly the precursors the search found rather than
 * every peptide at every charge in a range.
 *
 * <p>That makes this class a filter rather than a range. A named charge outside 1-6 is a <em>hard
 * error</em> that fails the whole {@code build-library} run: the models encode charge as a six-way
 * one-hot, so anything else would be predicted from no charge information at all, and FragCast
 * refuses rather than doing so. Those precursors are therefore left out of the list and counted
 * here, and {@link peptideptmformatting.PeptideSkipper} treats them as unsupported - the same way
 * MSBooster handles every other model's limits.
 */
public final class FragCastCharges {
    /** The charges FragCast's one-hot can represent. */
    public static final int LOWEST = 1;
    public static final int HIGHEST = 6;

    //pin files are read on several threads, so the tally is atomic
    private static final AtomicInteger unrepresentable = new AtomicInteger(0);
    private static final AtomicInteger reported = new AtomicInteger(0);

    private FragCastCharges() {
    }

    /**
     * Can FragCast predict this precursor? A charge it cannot represent is counted, so the run can
     * say how many PSMs it left unpredicted rather than simply producing fewer than expected.
     */
    public static boolean canPredict(int charge) {
        if (charge < LOWEST || charge > HIGHEST) {
            unrepresentable.incrementAndGet();
            return false;
        }
        return true;
    }

    /** As {@link #canPredict(int)}, for a charge that arrives as text. Unparseable is unusable. */
    public static boolean canPredict(String charge) {
        try {
            return canPredict(Integer.parseInt(charge.trim()));
        } catch (NumberFormatException e) {
            unrepresentable.incrementAndGet();
            return false;
        }
    }

    /**
     * How many PSMs (pin path) or precursors (peptide-list path) were left out because FragCast
     * cannot represent their charge. One count per {@link #canPredict} call, so it is a PSM count
     * when the charges came off pin rows and a precursor count when they came off a peptide list.
     */
    public static int unrepresentableCount() {
        return unrepresentable.get();
    }

    /**
     * Say once how many PSMs or precursors were left out. The generic "PSMs whose precursors are
     * not in the predictions" tally at the end of a run would otherwise be the only trace, and it
     * does not say why.
     *
     * @param unit what the tally counted: "PSMs" when {@link #canPredict} was asked once per pin
     *             row, "precursors" when it was asked once per peptide-list row
     */
    public static void reportSkipped(String unit) {
        int skipped = unrepresentableCount();
        if (skipped > 0 && reported.compareAndSet(0, 1)) {
            printInfo(skipped + " " + unit + " carry a precursor charge outside " + LOWEST + "-" + HIGHEST +
                    ", which FragCast's models cannot represent; they are left unpredicted");
        }
    }

    /** Forget the tally. Only needed between jobs in one JVM, i.e. in tests. */
    public static void reset() {
        unrepresentable.set(0);
        reported.set(0);
    }
}
