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

package allconstants;

import java.util.HashMap;

public class LowercaseModelMapper {
    private static HashMap<String, String> lowercaseToModel = new HashMap<>();
    public LowercaseModelMapper() {
        lowercaseToModel.put("", "");
        lowercaseToModel.put("dia-nn", "DIA-NN");
        lowercaseToModel.put("alphapeptdeep", "alphapeptdeep");

        for (String model : ModelCollections.KoinaModels) {
            lowercaseToModel.put(model.toLowerCase(), model);
        }
    }

    public HashMap<String, String> getLowercaseToModel() {
        return lowercaseToModel;
    }
}
