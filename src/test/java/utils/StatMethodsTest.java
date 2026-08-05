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
import static org.junit.jupiter.api.Assertions.assertFalse;

import java.util.concurrent.ConcurrentSkipListMap;

import org.junit.jupiter.api.Test;

/**
 * Covers {@link StatMethods#lookupInverse}, which reads the RT calibration
 * table {@code MzmlReader.irtToMinutes} — a predicted-RT to experimental-minutes
 * map — to produce {@code predRTrealUnits}.
 */
public class StatMethodsTest {

    private static ConcurrentSkipListMap<Double, Double> table() {
        ConcurrentSkipListMap<Double, Double> m = new ConcurrentSkipListMap<>();
        m.put(100.0, 50.0);
        m.put(146.205, 90.0);
        m.put(169.837, 99.2);
        return m;
    }

    @Test
    public void interpolatesBetweenBracketingKeys() {
        // Halfway between 146.205 and 169.837 is halfway between 90.0 and 99.2.
        double midKey = 0.5 * (146.205 + 169.837);
        assertEquals(0.5 * (90.0 + 99.2), StatMethods.lookupInverse(table(), midKey), 1e-9);
    }

    @Test
    public void exactKeyMatchReturnsThatKeysValue() {
        // floorEntry and ceilingEntry return the same entry here, so a naive
        // interpolation divides zero by zero and yields NaN. A NaN escapes both
        // the FeatureCalculator and PinWriter range gates, which then disagree
        // about whether the PSM was scored, and the pin writer throws.
        assertEquals(99.2, StatMethods.lookupInverse(table(), 169.837), 1e-9);
        assertEquals(50.0, StatMethods.lookupInverse(table(), 100.0), 1e-9);
        assertEquals(90.0, StatMethods.lookupInverse(table(), 146.205), 1e-9);
    }

    @Test
    public void exactKeyMatchIsNeverNaN() {
        for (double key : new double[]{100.0, 146.205, 169.837}) {
            assertFalse(Double.isNaN(StatMethods.lookupInverse(table(), key)),
                    "lookupInverse returned NaN for exact key " + key);
        }
    }

    @Test
    public void clampsOutsideTheTableRange() {
        assertEquals(50.0, StatMethods.lookupInverse(table(), 10.0), 1e-9);
        assertEquals(99.2, StatMethods.lookupInverse(table(), 500.0), 1e-9);
    }

    @Test
    public void singleEntryTableReturnsItsOnlyValue() {
        ConcurrentSkipListMap<Double, Double> m = new ConcurrentSkipListMap<>();
        m.put(120.0, 60.0);
        assertEquals(60.0, StatMethods.lookupInverse(m, 120.0), 1e-9);
        assertEquals(60.0, StatMethods.lookupInverse(m, 5.0), 1e-9);
        assertEquals(60.0, StatMethods.lookupInverse(m, 900.0), 1e-9);
    }
}
