package readers;

import readers.datareaders.ProteinHeaderParser;

public class ProteinHeaderParserTest {
    
    private static int passCount = 0;
    private static int failCount = 0;
    
    private static void assertEquals(Object expected, Object actual) {
        if ((expected == null && actual == null) || (expected != null && expected.equals(actual))) {
            passCount++;
        } else {
            failCount++;
            System.err.println("FAILED: Expected [" + expected + "] but got [" + actual + "]");
            StackTraceElement[] stack = Thread.currentThread().getStackTrace();
            System.err.println("  at " + stack[2]);
        }
    }
    
    private static void assertTrue(boolean condition) {
        if (condition) {
            passCount++;
        } else {
            failCount++;
            System.err.println("FAILED: Expected true but got false");
            StackTraceElement[] stack = Thread.currentThread().getStackTrace();
            System.err.println("  at " + stack[2]);
        }
    }
    
    private static void assertFalse(boolean condition) {
        if (!condition) {
            passCount++;
        } else {
            failCount++;
            System.err.println("FAILED: Expected false but got true");
            StackTraceElement[] stack = Thread.currentThread().getStackTrace();
            System.err.println("  at " + stack[2]);
        }
    }
    
    public static void testClassifyUniprot() {
        String header = "sp|P12345|PROT_HUMAN Protein description";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.UNIPROT, type);
    }
    
    public static void testClassifyNCBI() {
        String header = "NP_001234567.1 tumor protein p53 [Homo sapiens]";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.NCBI, type);
    }
    
    public static void testClassifyEnsembl() {
        String header = "ENSP00000269305 pep chromosome:GRCh38";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.ENSEMBL, type);
    }
    
    public static void testClassifyCptacEnsembl() {
        String header = "ENSP00000379334|ENST00000395929|ENSG00000141510 Gene";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.CPTAC_ENSEMBL, type);
    }
    
    public static void testClassifyNextprot() {
        String header = "nxp|NX_P04637|TP53 Tumor protein p53";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.NEXTPROT, type);
    }
    
    public static void testClassifyWithDecoyTag() {
        String header = "rev_sp|P12345|PROT_HUMAN Protein description";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.UNIPROT, type);
    }
    
    public static void testClassifyGeneric() {
        String header = "SomeRandomProtein12345 Description";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.GENERIC, type);
    }
    
    public static void testGetGeneNameUniprot() {
        String header = "sp|P12345|PROT_HUMAN Protein description GN=BRCA1 PE=1 SV=2";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.UNIPROT);
        assertEquals("BRCA1", geneName);
    }
    
    public static void testGetGeneNameNCBI() {
        String header = "NP_001234567.1 tumor protein p53 GN=TP53 [Homo sapiens]";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.NCBI);
        assertEquals("TP53", geneName);
    }
    
    public static void testGetGeneNameEnsembl() {
        String header = "ENSP00000269305 pep chromosome:GRCh38:17:43044295:43125483:1 gene:ENSG00000141510";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.ENSEMBL);
        assertEquals("ENSG00000141510", geneName);
    }
    
    public static void testGetGeneNameCptacEnsembl() {
        String header = "ENSP00000379334|ENST00000395929|ENSG00000141510.17|TP53";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.CPTAC_ENSEMBL);
        assertEquals("ENSG00000141510.17", geneName);
    }
    
    public static void testGetGeneNameNextprot() {
        String header = "nxp|NX_P04637|TP53|Tumor protein p53";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.NEXTPROT);
        assertEquals("TP53", geneName);
    }
    
    public static void testGetGeneNameTair() {
        String header = "AT1G01010 | Symbols: NAC001 | NAC domain protein | chr1:3760-5630";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.TAIR);
        assertEquals("NAC001", geneName);
    }
    
    public static void testGetGeneNameGeneric() {
        String header = "RandomProtein Description";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.GENERIC);
        assertEquals("", geneName);
    }
    
    public static void testGetGeneNameNoMatch() {
        String header = "sp|P12345|PROT_HUMAN Protein description without gene";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.UNIPROT);
        assertEquals("", geneName);
    }
    
    public static void testGetIDUniprot() {
        String header = "sp|P12345|PROT_HUMAN Protein description";
        String id = ProteinHeaderParser.getID(header, ProteinHeaderParser.DatabaseType.UNIPROT);
        assertEquals("P12345", id);
    }
    
    public static void testGetIDNCBI() {
        String header = "NP_001234567.1 tumor protein p53";
        String id = ProteinHeaderParser.getID(header, ProteinHeaderParser.DatabaseType.NCBI);
        assertEquals("NP_001234567.1", id);
    }
    
    public static void testGetIDEnsembl() {
        String header = "ENSP00000269305 pep chromosome:GRCh38";
        String id = ProteinHeaderParser.getID(header, ProteinHeaderParser.DatabaseType.ENSEMBL);
        assertEquals("ENSP00000269305", id);
    }
    
    public static void testGetIDNextprot() {
        String header = "nxp|NX_P04637|TP53 Tumor protein p53";
        String id = ProteinHeaderParser.getID(header, ProteinHeaderParser.DatabaseType.NEXTPROT);
        assertEquals("NX_P04637", id);
    }
    
    public static void testGetIDGeneric() {
        String header = "RandomProteinHeader";
        String id = ProteinHeaderParser.getID(header, ProteinHeaderParser.DatabaseType.GENERIC);
        assertEquals("RandomProteinHeader", id);
    }
    
    public static void testProcessHeaderUniprot() {
        String header = "sp|P12345|PROT_HUMAN Protein description GN=BRCA1 PE=1 SV=2";
        ProteinHeaderParser.ProteinRecord record = ProteinHeaderParser.processHeader(header, "rev_");
        
        assertEquals(header, record.originalHeader);
        assertEquals("sp|P12345|PROT_HUMAN", record.partHeader);
        assertEquals(ProteinHeaderParser.DatabaseType.UNIPROT, record.dbType);
        assertEquals("P12345", record.id);
        assertEquals("BRCA1", record.geneName);
        assertFalse(record.isDecoy);
    }
    
    public static void testProcessHeaderWithDecoy() {
        String header = "rev_sp|P12345|PROT_HUMAN Decoy protein GN=BRCA1";
        ProteinHeaderParser.ProteinRecord record = ProteinHeaderParser.processHeader(header, "rev_");
        
        assertEquals(ProteinHeaderParser.DatabaseType.UNIPROT, record.dbType);
        assertEquals("P12345", record.id);
        assertEquals("BRCA1", record.geneName);
        assertTrue(record.isDecoy);
    }
    
    public static void testProcessHeaderNCBI() {
        String header = "NP_000537.3 cellular tumor antigen p53 GN=TP53 [Homo sapiens]";
        ProteinHeaderParser.ProteinRecord record = ProteinHeaderParser.processHeader(header, "rev_");
        
        assertEquals(ProteinHeaderParser.DatabaseType.NCBI, record.dbType);
        assertEquals("NP_000537.3", record.id);
        assertEquals("TP53", record.geneName);
        assertFalse(record.isDecoy);
    }
    
    public static void testProcessHeaderEnsembl() {
        String header = "ENSP00000269305 pep chromosome:GRCh38:17:43044295:43125483:1 gene:ENSG00000141510";
        ProteinHeaderParser.ProteinRecord record = ProteinHeaderParser.processHeader(header, "rev_");
        
        assertEquals(ProteinHeaderParser.DatabaseType.ENSEMBL, record.dbType);
        assertEquals("ENSP00000269305", record.id);
        assertEquals("ENSG00000141510", record.geneName);
        assertFalse(record.isDecoy);
    }
    
    public static void testProcessHeaderNoSpaces() {
        String header = "sp|P12345|PROT_HUMAN";
        ProteinHeaderParser.ProteinRecord record = ProteinHeaderParser.processHeader(header, "rev_");
        
        assertEquals(header, record.partHeader);
    }
    
    public static void testMultipleGeneNamesUniprot() {
        String header = "sp|P12345|PROT_HUMAN Protein GN=GENE1 GN=GENE2 PE=1";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.UNIPROT);
        assertEquals("GENE1", geneName);
    }
    
    public static void testEnsemblGeneWithVersion() {
        String header = "ENSP00000379334 gene:ENSG00000141510.17 description";
        String geneName = ProteinHeaderParser.getGeneName(header, ProteinHeaderParser.DatabaseType.ENSEMBL);
        assertEquals("ENSG00000141510.17", geneName);
    }
    
    public static void testClassifyWithContaminantTag() {
        String header = "contam_sp|P12345|PROT_HUMAN Contaminant protein";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.UNIPROT, type);
    }
    
    public static void testClassifyWithDecoyAndContaminantTags() {
        String header = "rev_contam_sp|P12345|PROT_HUMAN Decoy contaminant";
        ProteinHeaderParser.DatabaseType type = ProteinHeaderParser.classify(header, "rev_");
        assertEquals(ProteinHeaderParser.DatabaseType.UNIPROT, type);
    }
    
    public static void main(String[] args) {
        System.out.println("Running ProteinHeaderParser tests...\n");
        
        testClassifyUniprot();
        testClassifyNCBI();
        testClassifyEnsembl();
        testClassifyCptacEnsembl();
        testClassifyNextprot();
        testClassifyWithDecoyTag();
        testClassifyGeneric();
        testClassifyWithContaminantTag();
        testClassifyWithDecoyAndContaminantTags();
        
        testGetGeneNameUniprot();
        testGetGeneNameNCBI();
        testGetGeneNameEnsembl();
        testGetGeneNameCptacEnsembl();
        testGetGeneNameNextprot();
        testGetGeneNameTair();
        testGetGeneNameGeneric();
        testGetGeneNameNoMatch();
        testMultipleGeneNamesUniprot();
        testEnsemblGeneWithVersion();
        
        testGetIDUniprot();
        testGetIDNCBI();
        testGetIDEnsembl();
        testGetIDNextprot();
        testGetIDGeneric();
        
        testProcessHeaderUniprot();
        testProcessHeaderWithDecoy();
        testProcessHeaderNCBI();
        testProcessHeaderEnsembl();
        testProcessHeaderNoSpaces();
        
        System.out.println("\n========================================");
        System.out.println("Test Results:");
        System.out.println("  Passed: " + passCount);
        System.out.println("  Failed: " + failCount);
        System.out.println("  Total:  " + (passCount + failCount));
        System.out.println("========================================");
        
        if (failCount == 0) {
            System.out.println("\nAll tests passed!");
            System.exit(0);
        } else {
            System.err.println("\nSome tests failed!");
            System.exit(1);
        }
    }
}

