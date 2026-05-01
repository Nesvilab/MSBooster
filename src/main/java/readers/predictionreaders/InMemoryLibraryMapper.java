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

package readers.predictionreaders;

import predictions.PredictionEntryHashMap;

/**
 * Adapter that exposes an existing {@link PredictionEntryHashMap} through the
 * {@link LibraryPredictionMapper} interface, used when predictions were produced
 * in-process (e.g. by FragPred) and there is no backing file to load from.
 */
public class InMemoryLibraryMapper implements LibraryPredictionMapper {
    private PredictionEntryHashMap allPreds;

    public InMemoryLibraryMapper(PredictionEntryHashMap allPreds) {
        this.allPreds = allPreds;
    }

    @Override
    public PredictionEntryHashMap getPreds() {
        return allPreds;
    }

    @Override
    public void setPreds(PredictionEntryHashMap preds) {
        this.allPreds = preds;
    }

    @Override
    public void clear() {
        if (allPreds != null) allPreds.clear();
    }
}
