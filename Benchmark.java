package texthack.app;

import texthack.algo.SuffixArrayIndex;
import texthack.bio.GenomeGenerator;
import texthack.bio.MotifHit;
import texthack.bio.MotifScanner;
import texthack.util.Lcp;
import texthack.util.Strs;

/**
 * Fair timing of the exact-matching algorithms plus suffix-array / LCP repeat analysis.
 *
 * Method: all algorithms get the SAME genome and the SAME prepared motifs. Genome loading and motif
 * validation are outside the timed region. The timed region is MotifScanner.run(), i.e. the algorithm's own
 * preprocessing (KMP failure table, Aho-Corasick automaton, hash set-up) plus the scan. A warm-up phase lets the
 * JIT compiler settle, then the search is repeated and the MEDIAN time is reported.
 */
public final class Benchmark {
    private Benchmark() {}

    /** Algorithms compared (suffix array is a different tool: index once, query many - see repeat analysis). */
    public static final String[] ALGOS = { "naive", "kmp", "rabin-karp", "aho-corasick", "z" };

    private static final class Timed {
        double medianMs;
        int runs;
        MotifHit[] hits;
    }

    private static Timed time(char[] text, MotifScanner.Job job, String algo) {
        Timed t = new Timed();
        long w0 = System.nanoTime(), now = w0;
        int warm = 0;
        MotifHit[] last = null;
        while (warm < 50 && (now - w0 < 100_000_000L || warm == 0)) {      // warm-up: up to ~100 ms
            last = MotifScanner.run(text, job, algo);
            warm++;
            now = System.nanoTime();
        }
        double avgMs = (now - w0) / 1e6 / warm;
        int reps = avgMs > 1000 ? 1 : avgMs > 200 ? 3 : 7;
        long[] ns = new long[reps];
        for (int i = 0; i < reps; i++) {
            long s = System.nanoTime();
            last = MotifScanner.run(text, job, algo);
            ns[i] = System.nanoTime() - s;
        }
        for (int i = 1; i < reps; i++) {                                    // insertion sort, tiny array
            long v = ns[i];
            int j = i - 1;
            while (j >= 0 && ns[j] > v) { ns[j + 1] = ns[j]; j--; }
            ns[j + 1] = v;
        }
        t.medianMs = ns[reps / 2] / 1e6;
        t.runs = reps;
        t.hits = last;
        return t;
    }

    /** Benchmarks all algorithms on one genome and prints the comparison table. Returns true if all agree. */
    public static boolean compare(char[] text, String[] motifs, boolean quietHeader) {
        MotifScanner.Job job = MotifScanner.prepare(motifs, false);
        int totalLen = 0;
        for (char[] p : job.pats) totalLen += p.length;
        if (!quietHeader) {
            System.out.println("Genome size : " + text.length + " bases");
            System.out.println("Motifs      : " + job.pats.length + " (total length " + totalLen + ") on the forward strand");
        }
        System.out.printf("%-14s | %11s | %16s | %s%n", "Algorithm", "Match Count", "Time (ms)", "Same positions as Naive?");
        System.out.println(Strs.repeat('-', 72));
        MotifHit[] ref = null;
        boolean allOk = true;
        for (String a : ALGOS) {
            Timed t = time(text, job, a);
            boolean ok = true;
            if (ref == null) ref = t.hits; else ok = MotifScanner.sameHits(ref, t.hits);
            allOk &= ok;
            System.out.printf("%-14s | %11d | %16.3f | %s%n", MotifScanner.displayName(a),
                    MotifScanner.totalMatches(t.hits), t.medianMs, ref == t.hits ? "(reference)" : ok ? "yes" : "NO  <-- MISMATCH");
        }
        return allOk;
    }

    /** Option 4: benchmark on the loaded genome, then a worst-case input for contrast. */
    public static void run(char[] text, String[] motifs) {
        boolean ok = compare(text, motifs, false);
        System.out.println();
        System.out.println(ok ? "Verification: every algorithm returned identical positions and counts."
                              : "Verification FAILED: algorithms disagree (this is a bug).");
        System.out.println("Timing: median of repeated runs after a warm-up; includes each algorithm's own preprocessing");
        System.out.println("        (e.g. building the Aho-Corasick automaton) but not file loading or input checks.");
        System.out.println("Note  : times vary between runs and computers, and no algorithm is fastest on every input.");
        System.out.println();
        worstCase();
    }

    /** Text AAAA...A with pattern AAAA...AC: naive re-compares almost the whole pattern at every shift. */
    public static void worstCase() {
        int n = 400_000, m = 600;
        char[] text = new char[n];
        for (int i = 0; i < n; i++) text[i] = 'A';
        String motif = Strs.repeat('A', m - 1) + "C";
        System.out.println("Worst-case input for naive: text = " + n + " x 'A', motif = " + (m - 1) + " x 'A' + 'C' (0 matches)");
        compare(text, new String[] { motif }, true);
    }

    /** Timing at increasing genome sizes (CLI: "scale"). One row per size, median ms per algorithm. */
    public static void scaling(int[] sizes, String[] motifs) {
        System.out.printf("%-12s", "Genome size");
        for (String a : ALGOS) System.out.printf(" | %12s", MotifScanner.displayName(a));
        System.out.printf(" | %s%n", "All agree?");
        System.out.println(Strs.repeat('-', 12 + ALGOS.length * 15 + 13));
        MotifScanner.Job job = MotifScanner.prepare(motifs, false);
        for (int n : sizes) {
            char[] g = GenomeGenerator.generate(n, 42, motifs, 50);
            System.out.printf("%-12d", n);
            MotifHit[] ref = null;
            boolean ok = true;
            for (String a : ALGOS) {
                Timed t = time(g, job, a);
                if (ref == null) ref = t.hits; else ok &= MotifScanner.sameHits(ref, t.hits);
                System.out.printf(" | %9.2f ms", t.medianMs);
            }
            System.out.printf(" | %s%n", ok ? "yes" : "NO");
        }
    }

    /** Option 5: suffix array + Kasai LCP, then the longest repeated region. */
    public static void repeatAnalysis(char[] text) {
        if (text.length > 5_000_000) System.out.println("Large genome: building the suffix array may take a while...");
        long t0 = System.nanoTime();
        SuffixArrayIndex idx = new SuffixArrayIndex(text);
        double saMs = (System.nanoTime() - t0) / 1e6;
        t0 = System.nanoTime();
        idx.lcp();
        double lcpMs = (System.nanoTime() - t0) / 1e6;
        t0 = System.nanoTime();
        int[] r = idx.longestRepeat();
        double scanMs = (System.nanoTime() - t0) / 1e6;
        System.out.printf("Suffix array built in %.1f ms | LCP array (Kasai) in %.1f ms | repeat search in %.1f ms%n", saMs, lcpMs, scanMs);
        if (r[2] == 0) {
            System.out.println("No repeated region found.");
            return;
        }
        String seq = new String(text, r[0], Math.min(r[2], 60)) + (r[2] > 60 ? "..." : "");
        System.out.println("Longest repeated region (stretches of N ignored):");
        System.out.println("  Length    : " + r[2] + " bases");
        System.out.println("  Sequence  : " + seq);
        System.out.println("  Positions : " + r[0] + " and " + r[1] + "   (0-based start; the two copies may overlap)");
    }
}
