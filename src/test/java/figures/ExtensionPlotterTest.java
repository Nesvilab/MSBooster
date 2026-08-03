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
