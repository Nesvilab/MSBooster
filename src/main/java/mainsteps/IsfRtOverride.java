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

package mainsteps;

import predictions.PredictionEntry;
import predictions.PredictionEntryHashMap;
import readers.datareaders.IsfAnnotation;

import java.util.concurrent.atomic.AtomicInteger;

import static utils.Print.printInfo;

/**
 * Finds the prediction whose RT a PSM that MSFragger marked as an in-source fragment (ISF) is
 * scored against, and counts what it did for one pin. Thread safe: PSMs of a pin are matched in
 * parallel.
 *
 * <p>An ISF ion is made in the source from a co-eluting parent ion, so it elutes with the parent.
 * Its expected RT is the parent's prediction, which exists because the parent is itself a pin row.
 * When it does not, or is 0 (no RT predicted), the PSM keeps its own prediction. In a mass-offset
 * search the parent's key also picks the RT calibration curve: its prediction lives on the parent's
 * offset group curve, not on the one of the ISF's own name.
 */
public final class IsfRtOverride {
    private final AtomicInteger matched = new AtomicInteger();
    private final AtomicInteger overridden = new AtomicInteger();

    /**
     * The parent prediction an ISF PSM is scored against, or null when it is missing or has no RT
     * (the PSM then keeps its own).
     */
    public PredictionEntry parentPrediction(IsfAnnotation isf, PredictionEntryHashMap allPreds) {
        matched.incrementAndGet();
        PredictionEntry parent = allPreds.get(isf.parentBaseCharge);
        if (parent == null || parent.RT == 0f) {
            return null;
        }
        overridden.incrementAndGet();
        return parent;
    }

    public int matched() {
        return matched.get();
    }

    public int overridden() {
        return overridden.get();
    }

    public int fallbacks() {
        return matched.get() - overridden.get();
    }

    /** One summary line per pin; nothing when the pin has no ISF annotations. */
    public void log(String pinName, int annotationsRead) {
        if (annotationsRead == 0) {
            return;
        }
        printInfo("In-source fragments in " + pinName + ": " + annotationsRead + " annotations read, "
                + matched() + " matched to pin rows, " + overridden() + " use the parent's predicted RT, "
                + fallbacks() + " fell back to their own predicted RT (parent prediction missing or 0)");
    }
}
