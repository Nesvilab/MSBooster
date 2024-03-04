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

package External;

import Features.Constants;
import Features.KoinaLibReader;
import Features.MainClass;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.Timer;
import java.util.TimerTask;
import java.util.concurrent.Callable;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.concurrent.atomic.AtomicLong;

public class KoinaTask implements Callable<Boolean> {
    private final String filename;
    private final String command;
    private final String property;
    private final String model;
    public int index;
    KoinaLibReader klr;
    AtomicLong waitTime;
    AtomicInteger numTimes;
    private int failedAttempts = 0;
    public boolean completed = false;

    public KoinaTask(String filename, String command, String property, String model, int index,
                     KoinaLibReader klr, AtomicLong waitTime, AtomicInteger numTimes) {
        this.filename = filename;
        this.command = command;
        this.property = property;
        this.model = model;
        this.index = index;
        this.klr = klr;
        this.waitTime = waitTime;
        this.numTimes = numTimes;
    }

    @Override
    public Boolean call() {
        StringBuilder koinaSb = new StringBuilder();
        try {
            //delay so we don't overwhelm server
            if (index < Constants.numThreads && failedAttempts == 0) {
                Thread.sleep(index * 100L);
            }

            long start = System.currentTimeMillis();
            ProcessBuilder builder = new ProcessBuilder(command.split(" "));
            builder.redirectErrorStream(true);
            Process process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));

            //Goal is to cut off 504 error early
            Timer timer = new Timer();
            TimerTask timerTask = new TimerTask() {
                @Override
                public void run() {
                    if (koinaSb.toString().isEmpty()) {
                        process.destroy();
                    }
                }
            };
            timer.schedule(timerTask, waitTime.get());

            String line = "";
            while ((line = reader.readLine()) != null) {
                koinaSb.append(line);
            }
            reader.close();
            process.waitFor();
            process.destroy();

            timer.cancel();
            KoinaModelCaller.parseKoinaOutput(filename, koinaSb.toString(),
                    property, model, klr);
            long timeDiff = System.currentTimeMillis() - start;
            long currentWaitTime = waitTime.get();
            waitTime.set(currentWaitTime + ((2 * timeDiff - currentWaitTime) /
                    (Math.min(Constants.numThreads, numTimes.getAndIncrement()))
            )); //local mean
            return true;
        } catch (Exception e) {
            String ending = koinaSb.substring(Math.max(0, koinaSb.toString().length() - 1000));
            failedAttempts++;

            if (failedAttempts == Constants.numKoinaAttempts) {
                System.out.println(command);
                System.out.println(filename + " had output that ended in: ");
                System.out.println(ending);
                System.out.println("Retried calling " + filename + " " + failedAttempts +
                        " times. This many be fixable by sending prediction requests to Koina at a slower rate " +
                        "by lowering the --numThreads parameter in the parameter file.");
                System.out.println("Exiting");
                System.exit(1);
            }

            return false;
        }
    }
}
