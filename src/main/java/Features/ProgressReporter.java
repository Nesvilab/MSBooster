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

package Features;

public class ProgressReporter {
    int linesRead = 1;
    int currentPercent;
    int iterations;
    long startTime;

    public ProgressReporter(int iterations) {
        this.iterations = iterations;
        currentPercent = Constants.loadingPercent;
//        startTime = System.nanoTime();
    }

    public int progress() {
        linesRead += 1;
        while (linesRead > iterations * currentPercent / 100 && currentPercent <= 100) {
//            long endTime = System.nanoTime();
//            System.out.print("..." + currentPercent + "% (" + (endTime - startTime) / 1000000000 + "sec)");
            if (currentPercent == 100) {
                System.out.println("..." + currentPercent + "%");
                currentPercent++;
                break;
            } else {
                System.out.print("..." + currentPercent + "%");
            }
//            startTime = System.nanoTime();
            currentPercent += Constants.loadingPercent;
            currentPercent = Math.min(currentPercent, 100);
        }
        return currentPercent;
    }
}
