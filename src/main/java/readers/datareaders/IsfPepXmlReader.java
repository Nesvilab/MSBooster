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

package readers.datareaders;

import javax.xml.stream.XMLInputFactory;
import javax.xml.stream.XMLStreamConstants;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamReader;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Collections;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.stream.Stream;

import static utils.Print.printError;
import static utils.Print.printInfo;

/**
 * Reads MSFragger's in-source fragment (ISF) annotations for a pin from the pepXML files MSFragger
 * wrote beside it: {@code <base>.pepXML} for DDA, {@code <base>_rank<N>.pepXML} for DIA and DDA+.
 * FragPipe names the pin {@code <base>.pin} in both cases (see its CmdMSBooster).
 *
 * <p>pepXML files can be gigabytes, so they are streamed with StAX and only the {@code isf_*} search
 * scores are kept. A file that cannot be read contributes nothing: ISF handling only refines the RT
 * features, so it must never stop MSBooster.
 */
public final class IsfPepXmlReader {
    static final String PARENT_PEPTIDE = "isf_parent_peptide";
    static final String PARENT_CHARGE = "isf_parent_charge";
    /** The rank {@link #locate} gives a file whose hits carry their real hit_rank. */
    static final int RANK_FROM_HIT = 0;

    private static final Pattern PIN_EXTENSION = Pattern.compile("(?i)\\.pin$");
    private static final int READ_BUFFER_BYTES = 1 << 20;

    private IsfPepXmlReader() {}

    /**
     * The ISF annotations of every PSM of this pin, keyed by {@link IsfAnnotation#key}; empty when
     * MSFragger wrote none. Unmodifiable.
     */
    public static Map<String, IsfAnnotation> readForPin(File pin) {
        Map<Path, Integer> sources;
        try {
            sources = locate(pin);
        } catch (IOException e) {
            printError("Could not list pepXML files next to " + pin + " (" + e.getMessage()
                    + "). Skipping in-source fragment annotations.");
            return Collections.emptyMap();
        }
        if (sources.isEmpty()) {
            printInfo("No pepXML next to " + pin.getName() + ", so no in-source fragment annotations");
            return Collections.emptyMap();
        }

        Map<String, IsfAnnotation> merged = new HashMap<>();
        sources.forEach((path, fixedRank) -> parse(path, fixedRank).forEach(merged::putIfAbsent));
        if (merged.isEmpty()) {
            printInfo("No in-source fragment annotations in the pepXML of " + pin.getName());
        } else {
            printInfo("Read " + merged.size() + " in-source fragment annotations for " + pin.getName()
                    + " from " + sources.size() + " pepXML file(s)");
        }
        return Collections.unmodifiableMap(merged);
    }

    /**
     * {@code <base>.pepXML} -> {@link #RANK_FROM_HIT} and every {@code <base>_rank<N>.pepXML} -> N,
     * in the pin's folder and in its listing order.
     */
    static Map<Path, Integer> locate(File pin) throws IOException {
        String base = PIN_EXTENSION.matcher(pin.getName()).replaceFirst("");
        Pattern pattern = Pattern.compile(Pattern.quote(base) + "(?:_rank([0-9]+))?\\.(?i:pepxml)");
        Path folder = pin.getAbsoluteFile().toPath().getParent();

        Map<Path, Integer> sources = new LinkedHashMap<>();
        try (Stream<Path> files = Files.list(folder)) {
            files.filter(Files::isRegularFile).forEach(path -> {
                Matcher matcher = pattern.matcher(path.getFileName().toString());
                if (matcher.matches()) {
                    sources.put(path, matcher.group(1) == null ? RANK_FROM_HIT : Integer.parseInt(matcher.group(1)));
                }
            });
        }
        return sources;
    }

    /**
     * The annotated hits of one pepXML, keyed by {@link IsfAnnotation#key}; empty when the file
     * cannot be read to its end.
     *
     * @param fixedRank the rank of every hit, or {@link #RANK_FROM_HIT} to use each hit's hit_rank
     */
    private static Map<String, IsfAnnotation> parse(Path pepXml, int fixedRank) {
        try (InputStream in = new BufferedInputStream(Files.newInputStream(pepXml), READ_BUFFER_BYTES)) {
            XMLInputFactory factory = XMLInputFactory.newInstance();
            factory.setProperty(XMLInputFactory.SUPPORT_DTD, false);
            factory.setProperty(XMLInputFactory.IS_SUPPORTING_EXTERNAL_ENTITIES, false);
            XMLStreamReader reader = factory.createXMLStreamReader(in);
            try {
                HitCollector collector = new HitCollector(fixedRank);
                while (reader.hasNext()) {
                    int event = reader.next();
                    if (event == XMLStreamConstants.START_ELEMENT) {
                        collector.start(reader);
                    } else if (event == XMLStreamConstants.END_ELEMENT) {
                        collector.end(reader.getLocalName());
                    }
                }
                if (collector.invalidHits > 0) {
                    printInfo("Skipped " + collector.invalidHits + " hits with incomplete in-source fragment "
                            + "annotations in " + pepXml.getFileName());
                }
                return collector.annotations;
            } finally {
                reader.close();
            }
        } catch (IOException | XMLStreamException | RuntimeException e) {
            printError("Could not read in-source fragment annotations from " + pepXml + " ("
                    + e.getMessage() + "). Ignoring this file.");
            return Collections.emptyMap();
        }
    }

    /** Streaming state: the current spectrum query and search hit, and what has been collected. */
    private static final class HitCollector {
        private final int fixedRank;
        private final Map<String, IsfAnnotation> annotations = new HashMap<>();
        private int invalidHits = 0;

        private String scan;
        /** The hit's own hit_rank; read only when the file name does not fix the rank. */
        private String hitRank;
        private String parentModifiedPeptide;
        private String parentCharge;

        HitCollector(int fixedRank) {
            this.fixedRank = fixedRank;
        }

        void start(XMLStreamReader reader) {
            switch (reader.getLocalName()) {
                case "spectrum_query":
                    scan = reader.getAttributeValue(null, "start_scan");
                    break;
                case "search_hit":
                    hitRank = fixedRank == RANK_FROM_HIT ? reader.getAttributeValue(null, "hit_rank") : null;
                    parentModifiedPeptide = null;
                    parentCharge = null;
                    break;
                case "search_score":
                    String name = reader.getAttributeValue(null, "name");
                    if (PARENT_PEPTIDE.equals(name)) {
                        parentModifiedPeptide = reader.getAttributeValue(null, "value");
                    } else if (PARENT_CHARGE.equals(name)) {
                        parentCharge = reader.getAttributeValue(null, "value");
                    }
                    break;
                default:
                    break;
            }
        }

        void end(String localName) {
            if (!"search_hit".equals(localName)) {
                return;
            }
            if (parentModifiedPeptide == null && parentCharge == null) {
                return; //not an ISF
            }
            try {
                int rank = fixedRank == RANK_FROM_HIT ? Integer.parseInt(hitRank) : fixedRank;
                IsfAnnotation annotation = new IsfAnnotation(parentModifiedPeptide, Integer.parseInt(parentCharge));
                annotations.put(IsfAnnotation.key(Integer.parseInt(scan), rank), annotation);
            } catch (IllegalArgumentException e) { //includes NumberFormatException
                invalidHits++;
            }
        }
    }
}
