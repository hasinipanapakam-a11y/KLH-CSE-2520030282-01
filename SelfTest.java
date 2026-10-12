package texthack.app;

import texthack.algo.SuffixArrayIndex;
import texthack.bio.*;
import texthack.util.*;
import java.io.*;

/** Dependency-free test runner: "Main test". Returns the number of failed checks. */
public final class SelfTest {
    private SelfTest() {}

    private static int passed, failed;

    private static void check(String name, boolean ok) {
        if (ok) passed++; else failed++;
        System.out.println((ok ? "  PASS  " : "  FAIL  ") + name);
    }

    private static String fmt(IntList l) {
        StringBuilder s = new StringBuilder();
        for (int i = 0; i < l.size(); i++) s.append(i > 0 ? "," : "").append(l.get(i));
        return s.toString();
    }

    /** Result of one motif for one algorithm as "p1,p2,..". */
    private static String find(String genome, String motif, String algo) {
        return fmt(MotifScanner.scan(genome.toCharArray(), new String[] { motif }, algo, false)[0].positions);
    }

    /** True if every algorithm returns exactly the expected position list. */
    private static boolean allAlgos(String genome, String motif, String expected) {
        for (String a : MotifScanner.ALGOS) if (!find(genome, motif, a).equals(expected)) {
            System.out.println("        " + a + " gave [" + find(genome, motif, a) + "] for " + motif + " in " + genome + ", expected [" + expected + "]");
            return false;
        }
        return true;
    }

    private static boolean throwsOn(String motif) {
        try { Dna.toMotif(motif); return false; } catch (IllegalArgumentException e) { return true; }
    }

    public static int run() {
        passed = 0; failed = 0;
        System.out.println("== Specification genome: TestGenome = ATGCGATAGCGATGCG (0-based positions) ==");
        String tg = "ATGCGATAGCGATGCG";
        check("ATG -> 0,11 (count 2), all algorithms", allAlgos(tg, "ATG", "0,11"));
        check("GCG -> 2,8,13 (count 3), all algorithms", allAlgos(tg, "GCG", "2,8,13"));
        check("TAG -> 6 (count 1), all algorithms", allAlgos(tg, "TAG", "6"));

        System.out.println("== Edge cases ==");
        check("motif that does not occur (CCCC) -> no matches", allAlgos(tg, "CCCC", ""));
        check("single-character motif A -> 0,5,7,11", allAlgos(tg, "A", "0,5,7,11"));
        check("overlapping matches: AAA in AAAAA -> 0,1,2", allAlgos("AAAAA", "AAA", "0,1,2"));
        check("overlapping matches: ATAT in ATATATAT -> 0,2,4", allAlgos("ATATATAT", "ATAT", "0,2,4"));
        check("motif at the very beginning (ATGCGA -> 0)", allAlgos(tg, "ATGCGA", "0"));
        check("motif at the very end (GATGCG -> 10, ends at the last base)", allAlgos(tg, "GATGCG", "10"));
        check("repeated motif at both ends (TGCG -> 1,12)", allAlgos(tg, "TGCG", "1,12"));
        check("motif equal to the whole genome -> 0", allAlgos(tg, tg, "0"));
        check("motif longer than genome -> no matches", allAlgos("ACGT", "ACGTA", ""));
        check("lower-case motif is accepted (atg -> 0,11)", allAlgos(tg, "atg", "0,11"));

        System.out.println("== Multiple motifs at once ==");
        MotifHit[] h = MotifScanner.scan(tg.toCharArray(), new String[] { "ATG", "GCG", "TAG", "CCCC" }, "aho-corasick", false);
        check("Aho-Corasick: 4 motifs in one pass", h.length == 4 && fmt(h[0].positions).equals("0,11") && fmt(h[1].positions).equals("2,8,13")
                && fmt(h[2].positions).equals("6") && h[3].positions.size() == 0);
        MotifScanner.Job dupJob = MotifScanner.prepare(new String[] { "ATG", "atg", "GCG" }, false);
        check("duplicate motifs are ignored", dupJob.pats.length == 2 && dupJob.duplicatesIgnored == 1);
        MotifHit[] over = MotifScanner.scan("AAAAA".toCharArray(), new String[] { "AAA", "AA" }, "aho-corasick", false);
        check("Aho-Corasick keeps overlapping matches of nested motifs", fmt(over[0].positions).equals("0,1,2") && fmt(over[1].positions).equals("0,1,2,3"));

        System.out.println("== Input validation ==");
        check("invalid base (ATX) rejected", throwsOn("ATX"));
        check("empty motif rejected", throwsOn("   "));
        check("digits rejected", throwsOn("A1G"));

        System.out.println("== Reverse complement (optional feature) ==");
        MotifHit[] rc = MotifScanner.scan("GTTTCC".toCharArray(), new String[] { "AAAC" }, "kmp", true);
        check("AAAC: reverse complement GTTT found on '-' strand at 0", rc.length == 2 && rc[1].strand == '-' && fmt(rc[1].positions).equals("0"));
        MotifHit[] pal = MotifScanner.scan("AGAATTCT".toCharArray(), new String[] { "GAATTC" }, "kmp", true);
        check("palindromic site GAATTC counted once (no double count)", pal.length == 1 && fmt(pal[0].positions).equals("1"));

        System.out.println("== FASTA loading and saving ==");
        try {
            File f = File.createTempFile("texthack", ".fasta");
            f.deleteOnExit();
            char[] src = GenomeGenerator.generate(1000, 5, new String[] { "ATG" }, 3);
            Fasta.write(f.getPath(), "TestGenome description", src);
            char[] back = Fasta.read(f.getPath());
            check("save then load returns the same 1000 bases", Dna.equal(src, back));
            check("header name read back as TestGenome", "TestGenome".equals(Fasta.firstHeader(f.getPath())));
            try (PrintWriter pw = new PrintWriter(new FileWriter(f))) {
                pw.println(">rec1"); pw.println("acgt"); pw.println("ACRY"); pw.println(">rec2"); pw.println("TTTT");
            }
            check("lower-case folded to upper-case; IUPAC codes (R, Y) become N", new String(Fasta.read(f.getPath())).equals("ACGTACNNTTTT"));
            try (PrintWriter pw = new PrintWriter(new FileWriter(f))) {
                pw.println(">rec1"); pw.println("ACGT"); pw.println(">rec2"); pw.println("TT");
            }
            check("two records are joined by a single N so motifs cannot span them", new String(Fasta.read(f.getPath())).equals("ACGTNTT"));
            boolean missing = false;
            try { Fasta.read("no_such_file.fasta"); } catch (IOException e) { missing = true; }
            check("missing file gives a clear IOException", missing);
        } catch (IOException e) {
            check("FASTA tests could run (" + e.getMessage() + ")", false);
        }

        System.out.println("== Randomised cross-check of all algorithms against Naive ==");
        Lcg r = new Lcg(2024);
        int bad = 0, trials = 3000;
        for (int t = 0; t < trials; t++) {
            int n = 1 + r.nextInt(150), k = 1 + r.nextInt(4);
            char[] tx = new char[n];
            for (int i = 0; i < n; i++) tx[i] = "ACGT".charAt(r.nextInt(k));
            int m = 1 + r.nextInt(7);
            char[] p = new char[m];
            for (int i = 0; i < m; i++) p[i] = "ACGT".charAt(r.nextInt(k));
            String[] one = { new String(p), "A" + new String(p) };
            for (boolean both : new boolean[] { false, true }) {
                MotifHit[] ref = MotifScanner.scan(tx, one, "naive", both);
                for (String a : MotifScanner.ALGOS) if (!MotifScanner.sameHits(ref, MotifScanner.scan(tx, one, a, both))) bad++;
            }
        }
        check(trials + " random texts x 2 motifs x 2 strand modes x 6 algorithms: " + bad + " mismatches", bad == 0);

        System.out.println("== Benchmark correctness ==");
        String[] motifs = { "ATG", "TATAAA", "GAATTC", "GGATCC", "AAGCTT", "GCCACCATGG" };
        char[] g = GenomeGenerator.generate(200_000, 42, motifs, 20);
        MotifHit[] ref = MotifScanner.scan(g, motifs, "naive", false);
        boolean same = true;
        for (String a : Benchmark.ALGOS) same &= MotifScanner.sameHits(ref, MotifScanner.scan(g, motifs, a, false));
        check("200,000-base genome: all benchmark algorithms return identical positions", same);
        check("planted motifs are all found (GCCACCATGG >= 10 forward copies)", ref[5].positions.size() >= 10);
        PrintStream real = System.out;
        boolean benchOk;
        try {
            System.setOut(new PrintStream(new ByteArrayOutputStream()));      // hide the table, keep only the verdict
            benchOk = Benchmark.compare(g, motifs, true);
        } finally { System.setOut(real); }
        check("Benchmark.compare reports that all algorithms agree", benchOk);

        System.out.println("== Suffix array, LCP and longest repeated region ==");
        char[] gat = "GATTACAGATTACA".toCharArray();
        int[] rep = new SuffixArrayIndex(gat).longestRepeat();
        check("GATTACAGATTACA: longest repeat GATTACA, length 7 at 0 and 7", rep[2] == 7 && rep[0] == 0 && rep[1] == 7);
        int[] repN = new SuffixArrayIndex("NNNNNNNNACGTAC".toCharArray()).longestRepeat();
        check("runs of N are ignored: NNNNNNNNACGTAC -> AC, length 2", repN[2] == 2);
        check("no repeat in ACGT -> length 0", new SuffixArrayIndex("ACGT".toCharArray()).longestRepeat()[2] == 0);
        boolean saOk = true;
        for (int t = 0; t < 300 && saOk; t++) {
            int n = 1 + r.nextInt(60), k = 1 + r.nextInt(4);
            char[] tx = new char[n];
            for (int i = 0; i < n; i++) tx[i] = "ACGT".charAt(r.nextInt(k));
            SuffixArrayIndex idx = new SuffixArrayIndex(tx);
            int[] sa = idx.suffixArray(), lcp = idx.lcp();
            for (int i = 1; i < n && saOk; i++) {
                int a = sa[i - 1], b = sa[i], l = 0;
                while (a + l < n && b + l < n && tx[a + l] == tx[b + l]) l++;
                boolean less = a + l == n || (b + l < n && tx[a + l] < tx[b + l]);       // suffix a < suffix b
                saOk = less && lcp[i] == l;
            }
        }
        check("300 random strings: suffix array sorted and Kasai LCP matches brute force", saOk);
        boolean sf = true;
        for (int t = 0; t < 300 && sf; t++) {
            int n = 1 + r.nextInt(80);
            char[] tx = new char[n];
            for (int i = 0; i < n; i++) tx[i] = "ACGT".charAt(r.nextInt(2));
            int[] best = new SuffixArrayIndex(tx).longestRepeat();
            int bl = 0;                                                    // brute-force longest repeat
            for (int a = 0; a < n; a++) for (int b = a + 1; b < n; b++) { int l = 0; while (b + l < n && tx[a + l] == tx[b + l]) l++; if (l > bl) bl = l; }
            sf = best[2] == bl;
        }
        check("300 random strings: longest-repeat length equals brute force", sf);

        System.out.println("== Increasing genome sizes (all algorithms must agree) ==");
        int[] sizes = { 1_000, 10_000, 100_000, 1_000_000 };
        for (int n : sizes) {
            char[] big = GenomeGenerator.generate(n, 42, motifs, 20);
            MotifHit[] base = MotifScanner.scan(big, motifs, "naive", true);
            boolean ok = true;
            for (String a : MotifScanner.ALGOS) ok &= MotifScanner.sameHits(base, MotifScanner.scan(big, motifs, a, true));
            check(n + " bases: all 6 algorithms agree (" + MotifScanner.totalMatches(base) + " matches)", ok);
        }

        System.out.println();
        System.out.println("Result: " + passed + " passed, " + failed + " failed");
        return failed;
    }
}
