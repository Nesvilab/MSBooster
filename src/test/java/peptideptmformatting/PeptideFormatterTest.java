/*
 * This file is part of MSBooster.
 *
 * MSBooster is free software: you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the Free
 * Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * MSBooster is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
 * for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MSBooster. If not, see <https://www.gnu.org/licenses/>.
 */

package peptideptmformatting;

import features.spectra.MassCalculator;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;

class PeptideFormatterTest {

    @Test
    void unsupportedCModificationKeepsCarbamidomethylatedProxyMass() {
        PeptideFormatter modified = new PeptideFormatter("AC[100.0000]DE", "2", "base");

        String modelPeptide = modified.getMs2pip();
        assertEquals("AC[UNIMOD:4]DE", modelPeptide);

        PeptideFormatter proxy = new PeptideFormatter(modelPeptide, "2", "ms2pip");
        assertEquals("AC[57.0215]DE", proxy.getBase());

        MassCalculator modifiedMc = new MassCalculator(modified.getBase(), "2");
        MassCalculator proxyMc = new MassCalculator(proxy.getBase(), "2");

        float b2Shift = modifiedMc.compareModMasses(proxyMc, 2, "b", 1, "");
        assertEquals(100.0 - PTMhandler.carbamidomethylationMass, b2Shift, 0.0001);

        float y2Shift = modifiedMc.compareModMasses(proxyMc, 2, "y", 1, "");
        assertEquals(0.0, y2Shift, 0.0001);
    }

    @Test
    void unsupportedNonCModificationIsStillRemovedFromModelInput() {
        PeptideFormatter modified = new PeptideFormatter("AM[100.0000]DE", "2", "base");

        assertEquals("AMDE", modified.getMs2pip());
    }

    @Test
    void predfullBareCProxyConvertsBackToCarbamidomethylatedBase() {
        PeptideFormatter modified = new PeptideFormatter("AC[100.0000]DE", "2", "base");

        String modelPeptide = modified.getPredfullKoina();
        assertEquals("ACDE", modelPeptide);

        PeptideFormatter proxy = new PeptideFormatter(modelPeptide, "2", "predfull");
        assertEquals("AC[57.0215]DE", proxy.getBase());

        MassCalculator modifiedMc = new MassCalculator(modified.getBase(), "2");
        MassCalculator proxyMc = new MassCalculator(proxy.getBase(), "2");

        float b2Shift = modifiedMc.compareModMasses(proxyMc, 2, "b", 1, "");
        assertEquals(100.0 - PTMhandler.carbamidomethylationMass, b2Shift, 0.0001);
    }
}
