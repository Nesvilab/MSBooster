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

import com.google.gson.JsonArray;
import com.google.gson.JsonElement;
import com.google.gson.JsonObject;
import com.google.gson.JsonParser;
import com.google.gson.JsonPrimitive;

import java.util.Arrays;
import java.util.Map;

/**
 * The precursors of one Koina request, in the order they are sent, which is the order of the rows Koina answers
 * with. KoinaModelCaller.parseKoinaOutput reads a response against them, so the response of any part of a request
 * can be read on its own.
 *
 * @param peptides the peptides, as written for the model
 * @param charges  the precursor charges, as written in the request; null for a model that is sent no charge
 */
public record KoinaBatch(String[] peptides, String[] charges) {
    public static final String PEPTIDE_SEQUENCES = "peptide_sequences";

    public KoinaBatch {
        if (peptides == null) {
            throw new IllegalArgumentException("A Koina batch needs its peptides");
        }
        if (charges != null && charges.length != peptides.length) {
            throw new IllegalArgumentException("A Koina batch of " + peptides.length + " peptides has "
                    + charges.length + " charges");
        }
    }

    public int size() {
        return peptides.length;
    }

    /** Precursors [from, to) of this batch. */
    public KoinaBatch slice(int from, int to) {
        return new KoinaBatch(Arrays.copyOfRange(peptides, from, to),
                charges == null ? null : Arrays.copyOfRange(charges, from, to));
    }

    /**
     * The request for precursors [from, to) of a request as JSONWriter writes it: every input cut to those rows,
     * with the first number of its shape set to how many there are. The values are written as the request wrote
     * them.
     *
     * @throws IllegalArgumentException if the request is not JSON with inputs of one row per precursor
     */
    public static String sliceRequest(String requestJson, int from, int to) {
        try {
            JsonObject request = JsonParser.parseString(requestJson).getAsJsonObject();
            JsonObject sliced = new JsonObject();
            for (Map.Entry<String, JsonElement> entry : request.entrySet()) {
                if (!entry.getKey().equals("inputs")) {
                    sliced.add(entry.getKey(), entry.getValue());
                    continue;
                }
                JsonArray inputs = new JsonArray();
                for (JsonElement element : entry.getValue().getAsJsonArray()) {
                    JsonObject input = element.getAsJsonObject().deepCopy();
                    JsonArray rows = input.getAsJsonArray("data");
                    if (to > rows.size()) {
                        throw new IllegalArgumentException("its input " + input.get("name") + " has " + rows.size()
                                + " rows, not " + to);
                    }
                    JsonArray cut = new JsonArray(to - from);
                    for (int i = from; i < to; i++) {
                        cut.add(rows.get(i));
                    }
                    input.add("data", cut);
                    JsonArray shape = input.getAsJsonArray("shape");
                    if (shape != null && shape.size() > 0) {
                        shape.set(0, new JsonPrimitive(to - from));
                    }
                    inputs.add(input);
                }
                sliced.add("inputs", inputs);
            }
            return sliced.toString();
        } catch (IllegalArgumentException e) {
            throw e;
        } catch (RuntimeException e) { //JsonParseException, or an element of another type than expected
            throw new IllegalArgumentException("Not a Koina request: " + e.getMessage(), e);
        }
    }

    /**
     * The precursors of a request as JSONWriter writes it: the rows of its peptide_sequences input, and of its
     * precursor_charges input if it has one.
     *
     * @throws IllegalArgumentException if the request is not JSON or has no peptide_sequences input
     */
    public static KoinaBatch fromRequest(String requestJson) {
        String[] peptides = null;
        String[] charges = null;
        try {
            JsonObject request = JsonParser.parseString(requestJson).getAsJsonObject();
            for (JsonElement element : request.getAsJsonArray("inputs")) {
                JsonObject input = element.getAsJsonObject();
                String name = input.get("name").getAsString();
                if (name.equals(PEPTIDE_SEQUENCES)) {
                    peptides = column(input);
                } else if (name.equals(KoinaModelSpec.Input.PRECURSOR_CHARGES.koinaName())) {
                    charges = column(input);
                }
            }
        } catch (RuntimeException e) { //JsonParseException, or an element of another type than expected
            throw new IllegalArgumentException("Not a Koina request: " + e.getMessage(), e);
        }
        if (peptides == null) {
            throw new IllegalArgumentException("A Koina request without " + PEPTIDE_SEQUENCES);
        }
        return new KoinaBatch(peptides, charges);
    }

    //the value of each row of an input shaped [n, 1], as written in the request
    private static String[] column(JsonObject input) {
        JsonArray rows = input.getAsJsonArray("data");
        String[] values = new String[rows.size()];
        for (int i = 0; i < values.length; i++) {
            JsonElement row = rows.get(i);
            values[i] = (row.isJsonArray() ? row.getAsJsonArray().get(0) : row).getAsString();
        }
        return values;
    }
}
