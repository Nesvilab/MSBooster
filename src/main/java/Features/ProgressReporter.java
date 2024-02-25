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

import java.util.concurrent.atomic.AtomicInteger;

public class ProgressReporter {
    private AtomicInteger linesRead = new AtomicInteger(1);
    private AtomicInteger currentPercent = new AtomicInteger(Constants.loadingPercent);
    private final int iterations;
    private boolean done = false;

    private static final Object lock = new Object(); // Create a lock object for synchronization

    public ProgressReporter(int iterations) {
        this.iterations = iterations;
    }

    public void progress() {
        synchronized (lock) {
            linesRead.incrementAndGet();
            while (linesRead.get() > (float) currentPercent.get() / 100f * (float) iterations &&
                    currentPercent.get() <= 100) {
                if (currentPercent.get() == 100) {
                    System.out.println("..." + currentPercent.get() + "%");
                    currentPercent.getAndIncrement();
                    done = true;
                    break;
                } else {
                    System.out.print("..." + currentPercent.get() + "%");
                }
                currentPercent.addAndGet(Constants.loadingPercent);
                currentPercent.set(Math.min(currentPercent.get(), 100));
            }
            if (linesRead.get() > iterations && !done) {
                System.out.println("..." + currentPercent.get() + "%");
            }
        }
    }
}
