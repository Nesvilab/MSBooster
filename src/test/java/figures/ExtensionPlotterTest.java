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

package figures;

import static figures.ExtensionPlotter.plot;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import allconstants.Constants;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.List;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.knowm.xchart.CategoryChart;
import org.knowm.xchart.CategoryChartBuilder;

class ExtensionPlotterTest {

    private static CategoryChart chartOf(List<Double> xData) {
        CategoryChart chart = new CategoryChartBuilder().width(600).height(400).build();
        chart.addSeries("targets", xData, List.of(1.0, 2.0));
        return chart;
    }

    // xchart opens the output file before it renders, so a chart that throws while painting used
    // to leave a 0 byte file behind. That empty figure reads as a real one and misdirects anyone
    // reading the stack trace, since the feature that actually failed looks like it succeeded.
    @Test
    void failedPlotLeavesNoPartialFile(@TempDir Path tmp) {
        Constants.plotExtension = "png";
        CategoryChart chart = chartOf(List.of(1.0, Double.POSITIVE_INFINITY));
        Path expected = tmp.resolve("broken.png");

        assertThrows(IllegalArgumentException.class, () -> plot(chart, tmp.resolve("broken").toString()));

        assertFalse(Files.exists(expected), "a failed render must not leave a file behind");
    }

    @Test
    void successfulPlotWritesFile(@TempDir Path tmp) throws Exception {
        Constants.plotExtension = "png";
        CategoryChart chart = chartOf(List.of(1.0, 2.0));
        Path expected = tmp.resolve("good.png");

        plot(chart, tmp.resolve("good").toString());

        assertTrue(Files.exists(expected), "a successful render must write its figure");
        assertTrue(Files.size(expected) > 0, "a successful render must not be empty");
    }
}
