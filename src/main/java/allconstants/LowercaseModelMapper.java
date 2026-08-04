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

package allconstants;

import java.util.HashMap;

public class LowercaseModelMapper {
    private static HashMap<String, String> lowercaseToModel = new HashMap<>();
    public LowercaseModelMapper() {
        lowercaseToModel.put("", "");
        lowercaseToModel.put("dia-nn", "DIA-NN");
        lowercaseToModel.put(FragCastModels.CONFORMER.toLowerCase(), FragCastModels.CONFORMER);
        lowercaseToModel.put(FragCastModels.FAST.toLowerCase(), FragCastModels.FAST);
        lowercaseToModel.put("alphapeptdeep", "alphapeptdeep");

        for (String model : ModelCollections.KoinaModels) {
            lowercaseToModel.put(model.toLowerCase(), model);
        }
    }

    public HashMap<String, String> getLowercaseToModel() {
        return lowercaseToModel;
    }
}
