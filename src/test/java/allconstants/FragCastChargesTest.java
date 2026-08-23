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

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

// FragCast expands only the list rows that omit a charge; a row that names one becomes exactly that
// precursor. MSBooster always names it, so this is a filter rather than a range - and the filter
// matters, because naming a charge FragCast cannot represent is not a dropped row but a hard error
// that fails the whole prediction.
public class FragCastChargesTest {
    @BeforeEach
    public void setUp() {
        FragCastCharges.reset();
    }

    @AfterEach
    public void tearDown() {
        FragCastCharges.reset();
    }

    @Test
    public void everyChargeTheOneHotCoversIsPredictable() {
        for (int charge = FragCastCharges.LOWEST; charge <= FragCastCharges.HIGHEST; charge++) {
            assertTrue(FragCastCharges.canPredict(charge), "charge " + charge + " was refused");
        }
        assertEquals(0, FragCastCharges.unrepresentableCount());
    }

    // Outside 1-6 the models have no charge encoding, and FragCast refuses the whole run rather
    // than predicting from none - so these must never reach the list.
    @Test
    public void chargesOutsideTheOneHotAreRefusedAndCounted() {
        for (int charge : new int[] {0, -1, 7, 8, 30}) {
            assertFalse(FragCastCharges.canPredict(charge), "charge " + charge + " was allowed");
        }
        assertEquals(5, FragCastCharges.unrepresentableCount());
    }

    @Test
    public void atextChargeIsReadAndAnUnreadableOneIsRefused() {
        assertTrue(FragCastCharges.canPredict("2"));
        assertTrue(FragCastCharges.canPredict(" 3 "), "surrounding whitespace made a charge unusable");
        assertFalse(FragCastCharges.canPredict("8"));
        assertFalse(FragCastCharges.canPredict(""));
        assertFalse(FragCastCharges.canPredict("2+"));
        assertEquals(3, FragCastCharges.unrepresentableCount());
    }

    // The count is what turns "fewer predictions than expected" into a stated reason.
    @Test
    public void theskippedCountIsReportedOnlyOnceAndOnlyWhenThereIsSomethingToReport() {
        FragCastCharges.reportSkipped();
        FragCastCharges.canPredict(9);
        assertEquals(1, FragCastCharges.unrepresentableCount());
        FragCastCharges.reportSkipped();
        FragCastCharges.reportSkipped();
    }
}
