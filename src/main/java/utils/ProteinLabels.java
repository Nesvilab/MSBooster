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

import java.util.LinkedHashSet;

/**
 * The one form in which a peptide's protein labels are handed to a predictor.
 *
 * <p>Both backends decide proteotypicity by asking a single question of this string - does it hold a
 * {@code ';'} - so every separator in it is a claim that the peptide sits in more than one protein.
 * That makes the separators load-bearing, and it makes a stray one a silent, whole-library error
 * rather than a cosmetic blemish: MSFragger's pin TERMINATES every entry, writing a peptide found in
 * one protein as {@code sp|P1|A_HUMAN;}, which read literally says "shared" for every peptide ever
 * searched. A library predicted from those labels comes back {@code Proteotypic = 0} on every row.
 *
 * <p>Hence one normaliser, called by every writer that fills that column, rather than a strip at each
 * site: the sources disagree with each other (MSFragger's pin terminates, its peptide list does not,
 * a hand-written list may do either), the predictors disagree about how forgiving to be, and a
 * peptide list may also come from a user's own file. Passing any of them through unexamined is what
 * put the defect in the library.
 *
 * <p>Only empty entries are dropped, and the ones that remain keep their order. A terminator is
 * exactly an empty last entry, so removing it needs no rule of its own, and the separators BETWEEN
 * real labels - the ones that carry the meaning - are left untouched.
 */
public final class ProteinLabels {
    private ProteinLabels() {}

    /**
     * The protein labels held in {@code fields}, as one {@code ';'}-joined list.
     *
     * <p>Several fields because a PIN may spread its labels across extra tab-separated columns its
     * header does not name; one field is the ordinary case. Duplicates collapse, which is what a
     * peptide listed once per column of a merged split-database search needs, and null or absent
     * fields are simply nothing to add - a pin with no protein column predicts an unlabelled library
     * rather than failing.
     */
    public static String normalize(String... fields) {
        if (fields == null) {
            return "";
        }
        final LinkedHashSet<String> labels = new LinkedHashSet<>();
        for (String field : fields) {
            if (field == null) {
                continue;
            }
            for (String label : field.split(";")) {
                final String trimmed = label.trim();
                if (!trimmed.isEmpty()) {
                    labels.add(trimmed);
                }
            }
        }
        return String.join(";", labels);
    }
}
