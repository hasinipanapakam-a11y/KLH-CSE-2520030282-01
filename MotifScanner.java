package texthack.bio;

import texthack.algo.*;
import texthack.util.Dna;
import texthack.util.IntList;

/**
 * Searches many motifs (optionally on both DNA strands) with a chosen algorithm.
 * Usage: prepare() validates the motifs once; run() performs only the search, so it can be timed fairly.
 * Positions are 0-based start offsets on the forward strand.
 */
public final class MotifScanner {
    /** Every selectable algorithm key (suffix-array is search-only, see Benchmark.ALGOS). */
    public static final String[] ALGOS = { "aho-corasick", "kmp", "rabin-karp", "z", "naive", "suffix-array" };

    private MotifScanner() {}

    /** Validated, de-duplicated search patterns (plus reverse complements when requested). */
    public static final class Job {
        public final char[][] pats;
        public final String[] names;
        public final char[] strands;
        public final int duplicatesIgnored;

        Job(char[][] pats, String[] names, char[] strands, int dup) {
            this.pats = pats; this.names = names; this.strands = strands; this.duplicatesIgnored = dup;
        }
    }

    public static String displayName(String algo) {
        switch (algo) {
            case "naive":        return "Naive";
            case "kmp":          return "KMP";
            case "z":            return "Z algorithm";
            case "rabin-karp":   return "Rabin-Karp";
            case "suffix-array": return "Suffix Array";
            case "aho-corasick": return "Aho-Corasick";
            default:             return algo;
        }
    }

    public static PatternMatcher matcher(String algo) {
        switch (algo) {
            case "naive":        return new NaiveMatcher();
            case "kmp":          return new KMPMatcher();
            case "z":            return new ZMatcher();
            case "rabin-karp":   return new RabinKarpMatcher();
            case "suffix-array": return new SuffixArrayMatcher();
            default: return null;   // aho-corasick is multi-pattern, handled in run()
        }
    }

    /** Validates motifs (A/C/G/T/N only), drops duplicates, adds reverse complements if asked. */
    public static Job prepare(String[] motifs, boolean bothStrands) {
        if (motifs == null || motifs.length == 0) throw new IllegalArgumentException("No motifs given");
        char[][] pats = new char[2 * motifs.length][];
        String[] names = new String[pats.length];
        char[] strands = new char[pats.length];
        int cnt = 0, dup = 0;
        for (String m : motifs) {
            char[] f = Dna.toMotif(m);
            boolean seen = false;
            for (int i = 0; i < cnt && !seen; i++) seen = strands[i] == '+' && Dna.equal(pats[i], f);
            if (seen) { dup++; continue; }
            String label = new String(f);
            pats[cnt] = f; names[cnt] = label; strands[cnt] = '+'; cnt++;
            if (bothStrands) {
                char[] rc = Dna.reverseComplement(f);
                if (!Dna.equal(rc, f)) {            // palindromic sites (e.g. GAATTC) would be double counted
                    pats[cnt] = rc; names[cnt] = label; strands[cnt] = '-'; cnt++;
                }
            }
        }
        char[][] p2 = new char[cnt][];
        String[] n2 = new String[cnt];
        char[] s2 = new char[cnt];
        for (int i = 0; i < cnt; i++) { p2[i] = pats[i]; n2[i] = names[i]; s2[i] = strands[i]; }
        return new Job(p2, n2, s2, dup);
    }

    /** Runs the search only (no validation). Returns one MotifHit per pattern, in job order. */
    public static MotifHit[] run(char[] genome, Job job, String algo) {
        int cnt = job.pats.length;
        IntList[] pos;
        if (algo.equals("aho-corasick")) {
            pos = new AhoCorasick(job.pats).search(genome);
        } else {
            PatternMatcher pm = matcher(algo);
            if (pm == null) throw new IllegalArgumentException("Unknown algorithm: " + algo);
            pos = new IntList[cnt];
            for (int i = 0; i < cnt; i++) pos[i] = pm.search(genome, job.pats[i]);
        }
        MotifHit[] out = new MotifHit[cnt];
        for (int i = 0; i < cnt; i++) out[i] = new MotifHit(job.names[i], job.strands[i], job.pats[i].length, pos[i]);
        return out;
    }

    public static MotifHit[] scan(char[] genome, String[] motifs, String algo, boolean bothStrands) {
        return run(genome, prepare(motifs, bothStrands), algo);
    }

    public static boolean sameHits(MotifHit[] a, MotifHit[] b) {
        if (a.length != b.length) return false;
        for (int i = 0; i < a.length; i++) if (!a[i].positions.sameAs(b[i].positions)) return false;
        return true;
    }

    public static long totalMatches(MotifHit[] h) {
        long t = 0;
        for (MotifHit x : h) t += x.positions.size();
        return t;
    }
}
