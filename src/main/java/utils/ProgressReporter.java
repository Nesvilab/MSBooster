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

import allconstants.Constants;

import java.util.concurrent.atomic.AtomicInteger;

public class ProgressReporter {
    private final AtomicInteger linesRead = new AtomicInteger(1);
    private final AtomicInteger currentPercent = new AtomicInteger(Constants.loadingPercent);
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
