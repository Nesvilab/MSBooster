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

package utils;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

// Every separator left in this string is a claim that the peptide sits in more than one protein -
// that is the whole of what a predictor reads it for - so these cases are about which separators
// survive, not about tidiness.
public class ProteinLabelsTest {
    // MSFragger's pin terminates every entry, so this is what a peptide in ONE protein looks like
    // coming out of a real pin. Left alone, it says "shared", and the predicted library comes back
    // Proteotypic = 0 on every row it ever wrote.
    @Test
    public void theTerminatorMsFraggerWritesIsNotASeparator() {
        assertEquals("sp|P12345|AAA_HUMAN", ProteinLabels.normalize("sp|P12345|AAA_HUMAN;"));
    }

    // The separators between real labels carry the meaning and must survive untouched, terminator
    // or no terminator - dropping one would call a shared peptide proteotypic instead.
    @Test
    public void separatorsBetweenRealLabelsSurvive() {
        assertEquals("sp|P1|A_HUMAN;sp|P2|B_HUMAN",
                ProteinLabels.normalize("sp|P1|A_HUMAN;sp|P2|B_HUMAN;"));
        assertEquals("sp|P1|A_HUMAN;sp|P2|B_HUMAN",
                ProteinLabels.normalize("sp|P1|A_HUMAN;sp|P2|B_HUMAN"));
    }

    // A pin may spread its labels across extra tab-separated columns its header does not name; the
    // peptide is just as shared as if they had arrived in one field.
    @Test
    public void labelsArrivingInSeveralFieldsBecomeOneList() {
        assertEquals("sp|P1|A_HUMAN;sp|P2|B_HUMAN",
                ProteinLabels.normalize("sp|P1|A_HUMAN;", "sp|P2|B_HUMAN;"));
    }

    // A merged split-database search can list the same protein once per chunk. Counting it twice
    // would make a proteotypic peptide read as shared.
    @Test
    public void aRepeatedLabelIsStillOneProtein() {
        assertEquals("sp|P1|A_HUMAN",
                ProteinLabels.normalize("sp|P1|A_HUMAN;", "sp|P1|A_HUMAN;"));
    }

    // Order is the predictor's choice of ProteinId: FragCast names the library after the first entry.
    @Test
    public void theFirstLabelStaysFirst() {
        assertEquals("sp|P2|B_HUMAN;sp|P1|A_HUMAN",
                ProteinLabels.normalize("sp|P2|B_HUMAN;sp|P1|A_HUMAN"));
    }

    // Nothing to label is a real answer: an unlabelled library beats a failed run, and a null column
    // must not reach a writer as the text "null".
    @Test
    public void nothingToLabelIsEmptyRatherThanText() {
        assertEquals("", ProteinLabels.normalize((String) null));
        assertEquals("", ProteinLabels.normalize((String[]) null));
        assertEquals("", ProteinLabels.normalize(""));
        assertEquals("", ProteinLabels.normalize(";"));
        assertEquals("", ProteinLabels.normalize("  ", null, ";;"));
    }

    // Blank entries anywhere in the list, not only the terminating one, are noise rather than
    // proteins - each would otherwise contribute a separator that says "shared".
    @Test
    public void blankEntriesAnywhereAreDropped() {
        assertEquals("sp|P1|A_HUMAN", ProteinLabels.normalize(";; sp|P1|A_HUMAN ;;"));
    }
}
