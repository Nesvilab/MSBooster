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

package koinaclasses;

import readers.predictionreaders.KoinaLibReader;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * Finds the precursors Koina refuses, so the others of their request are still predicted. Koina answers HTTP 400 to
 * a request with a precursor its model cannot take (KoinaTask.Result.REFUSED). Such a request is not sent again as
 * it is: it is cut in halves, each half is sent, and a refused half is cut again, until each precursor Koina refuses
 * is alone in a request. Those precursors are recorded on the KoinaLibReader ({@link KoinaLibReader#refuse}): they
 * have no prediction from the model, and transferKoinaPreds scores them as precursors the model cannot take.
 * <p>
 * Koina also answers HTTP 400 for a failure of its own, such as a model instance out of memory or reloading, which
 * the next request may not meet. An answer that names such a failure (CUDA, memory) never comes here: KoinaTask
 * counts it as a response not received, and the request is sent again whole. For a failure the answer does not name,
 * a request of one precursor is sent once more, {@link #CONFIRM_PAUSE_MS} after it was refused
 * ({@link KoinaTask#confirmation}), and its precursor counts as refused only when Koina refuses it again.
 * <p>
 * Koina also answers HTTP 400 when it cannot take any request of the model, for example for an input it does not
 * know. Cutting would then send 2n - 1 requests for each request of n precursors. So the first request refused
 * while nothing has been predicted is not cut at once: up to {@link #PROBES} of its precursors, spread over it, are
 * sent alone first, and the other refused requests wait. Once those precursors are back, the waiting requests are cut
 * if anything was predicted. If nothing was, the decision waits for the other requests still out, since one of them
 * may be predicted, and only when every request is back with nothing predicted does Koina refuse the model whatever
 * it is sent ({@link #refusesAll()}).
 * <p>
 * KoinaModelCaller calls it from the one thread that runs its requests, so it is not thread-safe.
 */
public final class KoinaBisection {
    /** The most precursors of a refused request sent alone to see whether Koina refuses every request. */
    public static final int PROBES = 8;
    /** The most characters of what Koina answered that a message quotes. */
    public static final int QUOTE_LENGTH = 300;
    /** How long a request of one precursor Koina refused waits before it is sent once more, in milliseconds. */
    public static final long CONFIRM_PAUSE_MS = 5000;
    //CONFIRM_PAUSE_MS, which tests shorten
    static volatile long confirmPauseMs = CONFIRM_PAUSE_MS;

    private enum State {
        /** no request was refused yet */
        NONE_REFUSED,
        /** a request was refused before anything was predicted, and some of its precursors are sent alone */
        PROBING,
        /** a refused request is cut in halves */
        CUTTING,
        /** Koina refuses every request of the model */
        REFUSES_ALL
    }

    /** A precursor Koina refused when it was sent alone: as the model was sent it, and what Koina answered. */
    public record Refusal(String peptide, String charge, String message) {}

    private final String model;
    private final KoinaLibReader klr;
    private State state = State.NONE_REFUSED;
    private boolean anyPredicted = false;
    //the requests other than probes that are not back yet: sent, waiting to be sent, or to be sent again
    private int pending;
    //the refused requests waiting for the probes to be back
    private final List<KoinaTask> held = new ArrayList<>();
    //the request the probes are taken from, and how many of them are not back yet
    private KoinaTask probed;
    private int probesSent = 0;
    private int probesOut = 0;
    private String firstRefusal;
    private final List<Refusal> isolated = new ArrayList<>();

    /**
     * @param requests the requests KoinaModelCaller sends first, one per request file
     */
    public KoinaBisection(String model, KoinaLibReader klr, int requests) {
        this.model = model;
        this.klr = klr;
        this.pending = requests;
    }

    /**
     * A request was predicted, or given up after its last attempt.
     *
     * @return the requests to send next
     */
    public List<KoinaTask> finished(KoinaTask task, boolean predicted) {
        if (predicted) {
            anyPredicted = true;
        }
        back(task);
        return decide();
    }

    /**
     * Koina refused a request (HTTP 400). A request of one precursor is sent once more, and recorded as refused when
     * Koina refuses it again; a longer one is cut in halves, or probed first (see the class).
     *
     * @return the requests to send next
     */
    public List<KoinaTask> refused(KoinaTask task) {
        if (firstRefusal == null) {
            firstRefusal = task.refusal();
        }
        if (task.size() == 1 && !task.isConfirmation()) {
            //it is not back yet: its confirmation stands for it
            return List.of(task.confirmation(confirmPauseMs));
        }
        back(task);
        List<KoinaTask> next = new ArrayList<>();
        if (task.size() == 1) {
            isolate(task);
        } else if (state == State.CUTTING || (state == State.NONE_REFUSED && anyPredicted)) {
            state = State.CUTTING;
            next.addAll(halves(task));
        } else if (state == State.NONE_REFUSED) {
            state = State.PROBING;
            List<KoinaTask> probes = probes(task);
            if (probes.size() < task.size()) {
                held.add(task);
            }
            next.addAll(probes);
        } else if (state == State.PROBING) {
            held.add(task);
        } //REFUSES_ALL: nothing more is sent

        next.addAll(decide());
        return next;
    }

    //a request is back for good: predicted, given up, or refused
    private void back(KoinaTask task) {
        if (task.isProbe()) {
            probesOut--;
        } else {
            pending--;
        }
    }

    //once the probes are back: the waiting requests are cut if anything was predicted. If nothing was, Koina
    //refuses every request once no other request is out that it may still predict
    private List<KoinaTask> decide() {
        if (state != State.PROBING || probesOut > 0) {
            return List.of();
        }
        if (!anyPredicted) {
            if (pending == 0) {
                state = State.REFUSES_ALL;
                held.clear();
            }
            return List.of();
        }
        state = State.CUTTING;
        List<KoinaTask> next = new ArrayList<>();
        for (KoinaTask task : held) {
            next.addAll(halves(task));
        }
        held.clear();
        return next;
    }

    private List<KoinaTask> halves(KoinaTask task) {
        int middle = task.size() / 2;
        pending += 2;
        return List.of(task.slice(0, middle, false), task.slice(middle, task.size(), false));
    }

    //up to PROBES precursors spread over the request, each alone
    private List<KoinaTask> probes(KoinaTask task) {
        int n = task.size();
        int count = Math.min(PROBES, n);
        List<KoinaTask> probes = new ArrayList<>(count);
        for (int i = 0; i < count; i++) {
            int precursor = (int) ((long) i * n / count);
            probes.add(task.slice(precursor, precursor + 1, true));
        }
        probed = task;
        probesSent = count;
        probesOut = count;
        return probes;
    }

    private void isolate(KoinaTask task) {
        KoinaBatch batch = task.batch();
        String charge = batch.charges() == null ? null : batch.charges()[0];
        if (klr.refuse(batch.peptides()[0], charge, task.refusal())) {
            isolated.add(new Refusal(batch.peptides()[0], charge, task.refusal()));
        }
    }

    /** Whether Koina refuses every request of the model: see the class. Nothing more is sent then. */
    public boolean refusesAll() {
        return state == State.REFUSES_ALL;
    }

    /** Why the model cannot predict the run, once {@link #refusesAll()}. */
    public String refusesAllMessage() {
        return "Koina model " + model + " refused every request MSBooster sent it, with HTTP 400 "
                + quote(firstRefusal) + ". It refused " + probed.label() + " (" + probed.size() + " precursors), and "
                + "each of the " + probesSent + " of its precursors sent alone, twice, so it cannot predict this "
                + "run's precursors. Check that the model takes this run's inputs (instrument, fragmentation, "
                + "collision energy)";
    }

    /** The precursors Koina refused when each was sent alone, in the order they were found. */
    public List<Refusal> isolated() {
        return Collections.unmodifiableList(isolated);
    }

    /** What Koina answered the first request it refused, as it can be quoted in a message; null if none. */
    public String firstRefusal() {
        return firstRefusal == null ? null : quote(firstRefusal);
    }

    /** What Koina answered, in one line, cut to {@value #QUOTE_LENGTH} characters. */
    public static String quote(String answer) {
        String line = answer == null ? "" : answer.replaceAll("\\s+", " ").trim();
        if (line.isEmpty()) {
            return "(no body)";
        }
        return line.length() > QUOTE_LENGTH ? line.substring(0, QUOTE_LENGTH) + "..." : line;
    }
}
