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

import static utils.Print.printError;
import static utils.Print.printInfo;

import Features.Constants;
import Features.KoinaLibReader;
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
    KoinaLibReader klr;
    AtomicLong waitTime;
    int failedAttempts = 0;
    public boolean completed = false;

    public KoinaTask(String filename, String command, String property, String model,
                     KoinaLibReader klr, AtomicLong waitTime) {
        this.filename = filename;
        this.command = command;
        this.property = property;
        this.model = model;
        this.klr = klr;
        this.waitTime = waitTime;
    }

    @Override
    public Boolean call() {
        StringBuilder koinaSb = new StringBuilder();
        Process process = null;
        try {
            long start = System.currentTimeMillis();
            ProcessBuilder builder = new ProcessBuilder(command.split(" "));
            builder.redirectErrorStream(true);
            process = builder.start();
            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream()));

            //Goal is to cut off 504 error early
            Timer timer = new Timer();
            Process finalProcess = process;
            TimerTask timerTask = new TimerTask() {
                @Override
                public void run() {
                    if (koinaSb.toString().isEmpty()) {
                        finalProcess.destroyForcibly();
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
            waitTime.set(currentWaitTime + ((3 * timeDiff - currentWaitTime) / 30));
            return true;
        } catch (Exception e) {
            try {
                String ending = koinaSb.substring(Math.max(0, koinaSb.toString().length() - 1000));
                process.destroyForcibly();

                if (failedAttempts == Constants.numKoinaAttempts) {
                    if (Constants.foundBest) {
                        printError(command);
                        printError(filename + " had output that ended in: ");
                        printError(ending);
                        printError("Retried calling " + filename + " " + failedAttempts + " times.");
                        printError("Exiting");
                        System.exit(1);
                    } else {
                        klr.failed = true;
                        return true;
                    }
                }

                return false;
            } catch (Exception e2) {
                e.printStackTrace();
                System.exit(1);
            }
            return false;
        }
    }
}
