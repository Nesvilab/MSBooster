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

package readers.datareaders;

import peptideptmformatting.PeptideFormatter;

/**
 * MSFragger's verdict that one PSM is an in-source fragment (ISF) of a co-eluting parent hit of the
 * same run, as written to its pepXML. Immutable.
 *
 * <p>An ISF ion is made in the source from the parent ion, so it elutes with the parent. The RT
 * MSBooster should expect for it is therefore the parent's predicted RT, which is looked up by
 * {@link #parentBaseCharge}.
 */
public final class IsfAnnotation {
    //MSFragger's pin Peptide column is prev.peptide+charge.next; the flanking residues are not
    //part of the key, so any placeholder works
    private static final String PIN_PREFIX = "-.";
    private static final String PIN_SUFFIX = ".-";

    /** The parent's modified peptide exactly as MSFragger writes it in the pin, minus flanks and charge. */
    public final String parentModifiedPeptide;
    public final int parentCharge;
    /** The prediction key MSBooster computes for the parent's own pin row. */
    public final String parentBaseCharge;

    /** Same "scan|rank" key {@code PepXMLDivider} uses to join pepXML hits to pin rows. */
    public static String key(int scanNum, int rank) {
        return scanNum + "|" + rank;
    }

    public IsfAnnotation(String parentModifiedPeptide, int parentCharge) {
        if (parentModifiedPeptide == null || parentModifiedPeptide.isEmpty()) {
            throw new IllegalArgumentException("ISF parent modified peptide must not be empty");
        }
        if (parentCharge < 1) {
            throw new IllegalArgumentException("ISF parent charge must be positive, got " + parentCharge);
        }
        this.parentModifiedPeptide = parentModifiedPeptide;
        this.parentCharge = parentCharge;
        //go through the pin formatter, exactly as for a pin row, so the key cannot drift from the
        //one the parent's prediction is stored under
        this.parentBaseCharge = new PeptideFormatter(
                PIN_PREFIX + parentModifiedPeptide + parentCharge + PIN_SUFFIX,
                String.valueOf(parentCharge), "pin").getBaseCharge();
    }
}
