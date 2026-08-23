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

package transferlearn.fragcast;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotNull;
import static org.junit.jupiter.api.Assertions.assertNull;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.io.OutputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import allconstants.Constants;
import allconstants.FragCastWeights;
import org.junit.jupiter.api.AfterEach;
import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;

/**
 * The bundle is the only thing standing between a user and three interchangeable-looking ONNX files.
 * These tests are mostly about what it refuses: the RT and IM models are the same architecture and
 * the same size, so a bundle that accepted a damaged or edited archive would predict the wrong
 * property without anything going wrong visibly.
 */
public class FragCastModelBundleTest {

    @TempDir
    Path tmp;

    private String outputDirectoryBefore;
    private String rtBefore;
    private String imBefore;
    private String specBefore;
    private String zipBefore;

    @BeforeEach
    public void setUp() {
        outputDirectoryBefore = Constants.outputDirectory;
        rtBefore = Constants.FragCastRtOnnx;
        imBefore = Constants.FragCastImOnnx;
        specBefore = Constants.FragCastSpecOnnx;
        zipBefore = Constants.FragCastModelZip;
        Constants.FragCastRtOnnx = "";
        Constants.FragCastImOnnx = "";
        Constants.FragCastSpecOnnx = "";
        Constants.FragCastModelZip = "";
    }

    @AfterEach
    public void tearDown() {
        Constants.outputDirectory = outputDirectoryBefore;
        Constants.FragCastRtOnnx = rtBefore;
        Constants.FragCastImOnnx = imBefore;
        Constants.FragCastSpecOnnx = specBefore;
        Constants.FragCastModelZip = zipBefore;
    }

    /** A stand-in for an ONNX file. Only the bytes matter to everything under test. */
    private File model(String name, String contents) throws IOException {
        File f = tmp.resolve(name).toFile();
        Files.write(f.toPath(), contents.getBytes(StandardCharsets.UTF_8));
        return f;
    }

    private FragCastWeights allThree() throws IOException {
        return FragCastWeights.of(
                model("rt.onnx", "retention time weights").getAbsolutePath(),
                model("im.onnx", "ion mobility weights").getAbsolutePath(),
                model("spec.onnx", "ms2 weights").getAbsolutePath());
    }

    // ---------------------------------------------------------------------------------------

    @Test
    public void roundTripsEveryModel() throws IOException {
        File zip = FragCastModelBundle.write(allThree(), tmp.resolve("bundle.zip").toFile());
        assertNotNull(zip);

        File dest = tmp.resolve("unpacked").toFile();
        FragCastWeights opened = FragCastModelBundle.open(zip, dest);

        // the entry names are the pretrained ones, so the unpacked directory can double as a
        // FRAGCAST_MODEL_DIR rather than only as a pile of files with a manifest beside them
        assertEquals(new File(dest, "FragCast-RT.onnx").getAbsolutePath(), opened.rtOnnx);
        assertEquals(new File(dest, "FragCast-IM.onnx").getAbsolutePath(), opened.imOnnx);
        assertEquals(new File(dest, "FragCast-Spec.onnx").getAbsolutePath(), opened.specOnnx);

        assertEquals("retention time weights", read(new File(opened.rtOnnx)));
        assertEquals("ion mobility weights", read(new File(opened.imOnnx)));
        assertEquals("ms2 weights", read(new File(opened.specOnnx)));
    }

    @Test
    public void carriesOnlyTheModelsItWasGiven() throws IOException {
        FragCastWeights partial = FragCastWeights.of(
                model("rt.onnx", "retention time weights").getAbsolutePath(), "", "");
        File zip = FragCastModelBundle.write(partial, tmp.resolve("rt-only.zip").toFile());

        FragCastWeights opened = FragCastModelBundle.open(zip, tmp.resolve("out").toFile());
        assertTrue(opened.imOnnx.isEmpty());
        assertTrue(opened.specOnnx.isEmpty());
        assertTrue(opened.rtOnnx.endsWith("FragCast-RT.onnx"));
    }

    @Test
    public void writesNothingWhenThereIsNothingToWrite() throws IOException {
        assertNull(FragCastModelBundle.write(FragCastWeights.base(), tmp.resolve("empty.zip").toFile()));
    }

    @Test
    public void refusesToBundleAModelThatIsNotThere() {
        FragCastWeights missing = FragCastWeights.of(tmp.resolve("gone.onnx").toString(), "", "");
        IOException e = assertThrows(IOException.class,
                () -> FragCastModelBundle.write(missing, tmp.resolve("x.zip").toFile()));
        assertTrue(e.getMessage().contains("no such file"), e.getMessage());
    }

    // --- what it refuses -------------------------------------------------------------------

    @Test
    public void refusesAnArchiveWithNoManifest() throws IOException {
        File zip = tmp.resolve("bare.zip").toFile();
        try (OutputStream out = Files.newOutputStream(zip.toPath());
             ZipOutputStream z = new ZipOutputStream(out)) {
            z.putNextEntry(new ZipEntry("FragCast-RT.onnx"));
            z.write("weights".getBytes(StandardCharsets.UTF_8));
            z.closeEntry();
        }
        IOException e = assertThrows(IOException.class,
                () -> FragCastModelBundle.open(zip, tmp.resolve("o1").toFile()));
        assertTrue(e.getMessage().contains(FragCastModelBundle.MANIFEST_NAME), e.getMessage());
    }

    @Test
    public void refusesAnArchiveWhoseContentsWereSwappedOut() throws IOException {
        File zip = FragCastModelBundle.write(allThree(), tmp.resolve("tampered.zip").toFile());
        // Rewrite one model to different content of exactly the same length, leaving the manifest's
        // checksum describing the file that used to be there. Equal length on purpose: it is what
        // gets past the size check and leaves the checksum as the only thing that can notice.
        rebuildWith(zip, "FragCast-IM.onnx", "someone elses wghts!");

        IOException e = assertThrows(IOException.class,
                () -> FragCastModelBundle.open(zip, tmp.resolve("o2").toFile()));
        assertTrue(e.getMessage().contains("checksum"), e.getMessage());
    }

    @Test
    public void refusesAnArchiveFromANewerMSBooster() throws IOException {
        File zip = FragCastModelBundle.write(allThree(), tmp.resolve("future.zip").toFile());
        rebuildManifest(zip, "bundleVersion=" + (FragCastModelBundle.BUNDLE_VERSION + 1) + "\ntasks=rt\n"
                + "rt.file=FragCast-RT.onnx\nrt.bytes=1\nrt.sha256=00\n");

        IOException e = assertThrows(IOException.class,
                () -> FragCastModelBundle.open(zip, tmp.resolve("o3").toFile()));
        assertTrue(e.getMessage().contains("version"), e.getMessage());
    }

    @Test
    public void refusesAnEntryThatWouldEscapeTheDestination() throws IOException {
        File zip = tmp.resolve("escape.zip").toFile();
        try (OutputStream out = Files.newOutputStream(zip.toPath());
             ZipOutputStream z = new ZipOutputStream(out)) {
            z.putNextEntry(new ZipEntry("../FragCast-RT.onnx"));
            z.write("weights".getBytes(StandardCharsets.UTF_8));
            z.closeEntry();
        }
        IOException e = assertThrows(IOException.class,
                () -> FragCastModelBundle.open(zip, tmp.resolve("o4").toFile()));
        assertTrue(e.getMessage().contains("plain file names"), e.getMessage());
    }

    @Test
    public void refusesAFileTheManifestDoesNotMention() throws IOException {
        // one model, so the extra entry stays under the entry-count cap and this test fails on the
        // manifest check it is named for rather than on the count
        FragCastWeights rtOnly = FragCastWeights.of(
                model("rt.onnx", "retention time weights").getAbsolutePath(), "", "");
        File zip = FragCastModelBundle.write(rtOnly, tmp.resolve("stowaway.zip").toFile());
        addEntry(zip, "notes.txt", "hello");

        IOException e = assertThrows(IOException.class,
                () -> FragCastModelBundle.open(zip, tmp.resolve("o5").toFile()));
        assertTrue(e.getMessage().contains("does not mention"), e.getMessage());
    }




    // --- wiring into a run -----------------------------------------------------------------

    @Test
    public void fillsInOnlyTheModelsTheParametersLeftEmpty() throws IOException {
        File zip = FragCastModelBundle.write(allThree(), tmp.resolve("apply.zip").toFile());
        Constants.outputDirectory = tmp.resolve("run").toString();
        Constants.FragCastModelZip = zip.getAbsolutePath();
        Constants.FragCastRtOnnx = "/somewhere/my-own-rt.onnx";

        FragCastModelBundle.applyToConstants();

        assertEquals("/somewhere/my-own-rt.onnx", Constants.FragCastRtOnnx);
        assertTrue(Constants.FragCastImOnnx.endsWith("FragCast-IM.onnx"), Constants.FragCastImOnnx);
        assertTrue(Constants.FragCastSpecOnnx.endsWith("FragCast-Spec.onnx"), Constants.FragCastSpecOnnx);
    }

    @Test
    public void doesNothingWhenNoBundleWasNamed() throws IOException {
        Constants.FragCastModelZip = "";
        FragCastModelBundle.applyToConstants();
        assertTrue(Constants.FragCastRtOnnx.isEmpty());

        // a params file that leaves the parameter at "null" means the same thing as leaving it out
        Constants.FragCastModelZip = "null";
        FragCastModelBundle.applyToConstants();
        assertTrue(Constants.FragCastRtOnnx.isEmpty());
    }

    @Test
    public void unpacksToOneStableDirectoryPerRun() throws IOException {
        // the paths are hashed into the prediction file name, so two calls in one run must agree or
        // every rerun re-predicts a library it already had
        Constants.outputDirectory = tmp.resolve("run").toString();
        File zip = tmp.resolve("x.zip").toFile();
        assertEquals(FragCastModelBundle.unpackDirectory(zip), FragCastModelBundle.unpackDirectory(zip));
        assertEquals(new File(Constants.outputDirectory, "fragcast_model"),
                FragCastModelBundle.unpackDirectory(zip));
    }

    // ---------------------------------------------------------------------------------------

    private static String read(File f) throws IOException {
        return new String(Files.readAllBytes(f.toPath()), StandardCharsets.UTF_8);
    }

    private String manifestText(File zip) throws IOException {
        try (java.util.zip.ZipFile z = new java.util.zip.ZipFile(zip)) {
            ZipEntry e = z.getEntry(FragCastModelBundle.MANIFEST_NAME);
            java.io.ByteArrayOutputStream buf = new java.io.ByteArrayOutputStream();
            try (java.io.InputStream in = z.getInputStream(e)) {
                byte[] chunk = new byte[4096];
                int n;
                while ((n = in.read(chunk)) > 0) {
                    buf.write(chunk, 0, n);
                }
            }
            return new String(buf.toByteArray(), StandardCharsets.UTF_8);
        }
    }

    private String manifestValue(File zip, String key) throws IOException {
        for (String line : manifestText(zip).split("\n")) {
            if (line.startsWith(key + "=")) {
                return line.substring(key.length() + 1).trim();
            }
        }
        return null;
    }

    /** Rewrite the archive with one entry's bytes replaced, leaving everything else alone. */
    private void rebuildWith(File zip, String entryName, String contents) throws IOException {
        rebuild(zip, entryName, contents.getBytes(StandardCharsets.UTF_8), false);
    }

    private void rebuildManifest(File zip, String manifest) throws IOException {
        rebuild(zip, FragCastModelBundle.MANIFEST_NAME, manifest.getBytes(StandardCharsets.UTF_8), false);
    }

    private void addEntry(File zip, String entryName, String contents) throws IOException {
        rebuild(zip, entryName, contents.getBytes(StandardCharsets.UTF_8), true);
    }

    private void rebuild(File zip, String entryName, byte[] contents, boolean append) throws IOException {
        java.util.LinkedHashMap<String, byte[]> entries = new java.util.LinkedHashMap<>();
        try (java.util.zip.ZipFile z = new java.util.zip.ZipFile(zip)) {
            java.util.Enumeration<? extends ZipEntry> it = z.entries();
            while (it.hasMoreElements()) {
                ZipEntry e = it.nextElement();
                java.io.ByteArrayOutputStream buf = new java.io.ByteArrayOutputStream();
                try (java.io.InputStream in = z.getInputStream(e)) {
                    byte[] chunk = new byte[4096];
                    int n;
                    while ((n = in.read(chunk)) > 0) {
                        buf.write(chunk, 0, n);
                    }
                }
                entries.put(e.getName(), buf.toByteArray());
            }
        }
        if (append || entries.containsKey(entryName)) {
            entries.put(entryName, contents);
        }
        try (OutputStream out = Files.newOutputStream(zip.toPath());
             ZipOutputStream z = new ZipOutputStream(out)) {
            for (java.util.Map.Entry<String, byte[]> e : entries.entrySet()) {
                z.putNextEntry(new ZipEntry(e.getKey()));
                z.write(e.getValue());
                z.closeEntry();
            }
        }
    }
}
