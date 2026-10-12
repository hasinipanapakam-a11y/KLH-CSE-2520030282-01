package texthack.app;

import texthack.bio.*;
import texthack.util.Dna;
import texthack.util.IntList;
import texthack.util.Strs;
import java.io.*;

/** Menu-driven front end of the DNA Sequence Analysis Tool. */
public final class Main {
    private static final BufferedReader IN = new BufferedReader(new InputStreamReader(System.in));
    private static final String[] DEFAULT_MOTIFS = { "ATG", "TATAAA", "GAATTC", "GGATCC", "AAGCTT", "GCCACCATGG" };
    private static final int PLANT_COPIES = 50;
    private static final int MIN_GENOME = 1_000, MAX_GENOME = 50_000_000;
    private static final int SHOW_POSITIONS = 10;

    private static char[] genome;
    private static String genomeName = "(none)";
    private static String[] lastMotifs = DEFAULT_MOTIFS;

    public static void main(String[] args) throws IOException {
        if (args.length > 0) {
            switch (args[0]) {
                case "test":  System.exit(SelfTest.run() == 0 ? 0 : 1); return;
                case "demo":  demo(); return;
                case "scale": scale(args); return;
                case "gen":
                    if (args.length < 3) { System.out.println("Usage: gen <bases> <file.fasta>"); return; }
                    char[] g = GenomeGenerator.generate(Integer.parseInt(args[1]), 42, DEFAULT_MOTIFS, PLANT_COPIES);
                    Fasta.write(args[2], "synthetic_genome", g);
                    System.out.println("Wrote " + g.length + " bases to " + args[2]);
                    return;
                default:
                    System.out.println("Unknown argument '" + args[0] + "'. Use: test | demo | scale | gen <bases> <file>");
                    return;
            }
        }
        while (true) {
            printMenu();
            String c = ask("Enter choice (1-7)", "");
            if (c == null) return;                       // input closed
            try {
                switch (c) {
                    case "1": loadGenome(); break;
                    case "2": generateGenome(); break;
                    case "3": if (hasGenome()) searchMotifs(); break;
                    case "4": if (hasGenome()) benchmark(); break;
                    case "5": if (hasGenome()) Benchmark.repeatAnalysis(genome); break;
                    case "6": if (hasGenome()) saveGenome(); break;
                    case "7": System.out.println("Goodbye."); return;
                    default:  System.out.println("Please enter a number from 1 to 7.");
                }
            } catch (IllegalArgumentException | IOException e) {
                System.out.println("Error: " + e.getMessage());
            }
        }
    }

    // ---------------------------------------------------------------- menu & input helpers

    private static void printMenu() {
        System.out.println();
        System.out.println("==============================");
        System.out.println("   DNA SEQUENCE ANALYSIS TOOL");
        System.out.println("==============================");
        System.out.println("Current genome: " + (genome == null ? "none loaded" : genomeName + " (" + genome.length + " bases)"));
        System.out.println();
        System.out.println("1. Load DNA Genome from FASTA");
        System.out.println("2. Generate Synthetic Genome");
        System.out.println("3. Search DNA Motifs");
        System.out.println("4. Compare Algorithms (Benchmark)");
        System.out.println("5. Find Longest Repeated Region");
        System.out.println("6. Save Genome as FASTA");
        System.out.println("7. Exit");
        System.out.println();
    }

    /** Prompts; returns the default when the user just presses Enter, or null when input has ended. */
    private static String ask(String q, String def) throws IOException {
        System.out.print(q + (def.isEmpty() ? "" : " [" + def + "]") + ": ");
        String s = IN.readLine();
        if (s == null) return null;
        s = s.trim();
        return s.isEmpty() ? def : s;
    }

    private static boolean hasGenome() {
        if (genome == null) System.out.println("No genome loaded yet. Choose option 1 (load) or 2 (generate) first.");
        return genome != null;
    }

    private static String stripQuotes(String s) {
        if (s.length() >= 2 && s.charAt(0) == '"' && s.charAt(s.length() - 1) == '"') return s.substring(1, s.length() - 1);
        return s;
    }

    private static String join(String[] a) {
        StringBuilder sb = new StringBuilder();
        for (int i = 0; i < a.length; i++) sb.append(i > 0 ? ", " : "").append(a[i]);
        return sb.toString();
    }

    // ---------------------------------------------------------------- 1, 2, 6: genome handling

    private static void loadGenome() throws IOException {
        String path = ask("FASTA file path", "sample/sample_genome.fasta");
        if (path == null) return;
        path = stripQuotes(path);
        char[] g = Fasta.read(path);
        String header = Fasta.firstHeader(path);
        genome = g;
        genomeName = header != null ? header : new File(path).getName();
        describeGenome();
    }

    private static void generateGenome() throws IOException {
        String s = ask("Genome length in bases", "1000000");
        if (s == null) return;
        String digits = s.replace(",", "").replace("_", "");
        if (!digits.matches("\\d{1,9}") || Long.parseLong(digits) < MIN_GENOME || Long.parseLong(digits) > MAX_GENOME)
            throw new IllegalArgumentException("Enter a whole number between " + MIN_GENOME + " and " + MAX_GENOME + ".");
        genome = GenomeGenerator.generate(Integer.parseInt(digits), 42, DEFAULT_MOTIFS, PLANT_COPIES);
        genomeName = "synthetic";
        describeGenome();
        System.out.println("Planted motifs (" + PLANT_COPIES + " copies each, on both strands): " + join(DEFAULT_MOTIFS));
    }

    private static void describeGenome() {
        long[] c = new long[5];
        for (char ch : genome) c[Dna.code(ch)]++;
        long acgt = c[0] + c[1] + c[2] + c[3];
        System.out.printf("Loaded '%s': %d bases | A=%d C=%d G=%d T=%d N=%d | GC content %.2f%%%n", genomeName, genome.length,
                c[0], c[1], c[2], c[3], c[4], acgt == 0 ? 0.0 : 100.0 * (c[1] + c[2]) / acgt);
    }

    private static void saveGenome() throws IOException {
        String path = ask("Output FASTA file", "genome.fasta");
        if (path == null) return;
        path = stripQuotes(path);
        Fasta.write(path, genomeName, genome);
        System.out.println("Saved " + genome.length + " bases to " + new File(path).getAbsolutePath());
    }

    // ---------------------------------------------------------------- 3: motif search

    private static void searchMotifs() throws IOException {
        System.out.println("Enter DNA motifs separated by commas (letters A, C, G, T, N). Press Enter for: " + join(lastMotifs));
        System.out.println("Optional extras at the end of the line: --rc (also search the reverse-complement strand), --csv [file] (export all positions)");
        String line = ask("Motifs", "");
        if (line == null) return;

        int dash = line.indexOf("--");
        String motifPart = dash < 0 ? line : line.substring(0, dash);
        boolean rc = false, csv = false;
        String csvPath = "motif_matches.csv";
        if (dash >= 0) {
            String[] tok = line.substring(dash).trim().split("\\s+");
            for (int i = 0; i < tok.length; i++) {
                if (tok[i].equals("--rc")) rc = true;
                else if (tok[i].equals("--csv")) {
                    csv = true;
                    if (i + 1 < tok.length && !tok[i + 1].startsWith("--")) csvPath = stripQuotes(tok[++i]);
                } else throw new IllegalArgumentException("Unknown option '" + tok[i] + "'. Valid options: --rc, --csv [file]");
            }
        }
        String[] motifs = lastMotifs;
        if (!motifPart.trim().isEmpty()) {
            String[] raw = Strs.splitComma(motifPart.replace(' ', ',').replace('\t', ','));
            int k = 0;
            for (String r : raw) if (!r.isEmpty()) k++;
            if (k == 0) throw new IllegalArgumentException("No motifs found in the input.");
            motifs = new String[k];
            k = 0;
            for (String r : raw) if (!r.isEmpty()) motifs[k++] = r;
        }
        MotifScanner.Job job = MotifScanner.prepare(motifs, rc);       // validates; throws a clear message

        System.out.println("Algorithms: 1) Aho-Corasick   2) KMP   3) Rabin-Karp   4) Z algorithm   5) Naive   6) Suffix Array");
        String algo = null;
        while (algo == null) {
            String a = ask("Algorithm", "1");
            if (a == null) return;
            if (a.matches("[1-6]")) algo = MotifScanner.ALGOS[Integer.parseInt(a) - 1];
            else System.out.println("Please enter a number from 1 to 6.");
        }

        long t0 = System.nanoTime();
        MotifHit[] hits = MotifScanner.run(genome, job, algo);
        double ms = (System.nanoTime() - t0) / 1e6;
        int plus = 0;
        for (char st : job.strands) if (st == '+') plus++;
        lastMotifs = new String[plus];                                 // cleaned list: upper-case, no duplicates
        plus = 0;
        for (int i = 0; i < job.names.length; i++) if (job.strands[i] == '+') lastMotifs[plus++] = job.names[i];
        printResults(hits, job, algo, ms, rc);
        if (csv) writeCsv(csvPath, hits);
    }

    private static void printResults(MotifHit[] hits, MotifScanner.Job job, String algo, double ms, boolean rc) {
        System.out.println();
        System.out.println("Genome    : " + genomeName + " (" + genome.length + " bases)");
        System.out.println("Algorithm : " + MotifScanner.displayName(algo));
        System.out.printf("Time      : %.3f ms (single run; use option 4 for repeated, more reliable timings)%n", ms);
        System.out.println("Positions : 0-based start index (the first base of the genome is position 0)");
        if (job.duplicatesIgnored > 0) System.out.println("Note      : " + job.duplicatesIgnored + " duplicate motif(s) ignored");
        System.out.println();
        System.out.printf("%-24s %s%-7s %s%n", "Motif", rc ? "Strand  " : "", "Count", "Positions");
        System.out.println(Strs.repeat('-', rc ? 78 : 70));
        for (MotifHit h : hits) {
            String label = h.motif.length() > 24 ? h.motif.substring(0, 21) + "..." : h.motif;
            System.out.printf("%-24s %s%-7d %s%n", label, rc ? String.format("%-8c", h.strand) : "", h.positions.size(), positions(h.positions));
        }
        boolean cut = false;
        for (MotifHit h : hits) cut |= h.positions.size() > SHOW_POSITIONS;
        if (cut) System.out.println("Only the first " + SHOW_POSITIONS + " positions per motif are shown; add --csv to export every position.");
        if (rc) System.out.println("Strand '-': the reverse complement of the motif was found; positions use forward-strand coordinates.");
    }

    private static String positions(IntList p) {
        if (p.size() == 0) return "not found";
        StringBuilder sb = new StringBuilder();
        int show = Math.min(SHOW_POSITIONS, p.size());
        for (int i = 0; i < show; i++) sb.append(i > 0 ? ", " : "").append(p.get(i));
        if (p.size() > show) sb.append(", ... (+").append(p.size() - show).append(" more)");
        return sb.toString();
    }

    private static void writeCsv(String path, MotifHit[] hits) throws IOException {
        long rows = 0;
        try (PrintWriter pw = new PrintWriter(new BufferedWriter(new FileWriter(path)))) {
            pw.println("motif,strand,start_0based,end_exclusive");
            for (MotifHit h : hits) {
                IntList p = h.positions;
                for (int i = 0; i < p.size(); i++) { pw.println(h.motif + "," + h.strand + "," + p.get(i) + "," + (p.get(i) + h.length)); rows++; }
            }
        }
        System.out.println("Exported " + rows + " match(es) to " + new File(path).getAbsolutePath());
    }

    // ---------------------------------------------------------------- 4: benchmark

    private static void benchmark() {
        System.out.println("Benchmarking with motifs: " + join(lastMotifs) + "  (change them via option 3)");
        System.out.println("Running warm-up and repeated timings, please wait...");
        System.out.println();
        Benchmark.run(genome, lastMotifs);
    }

    // ---------------------------------------------------------------- command-line extras

    private static void demo() {
        System.out.println("Generating a 1,000,000-base synthetic genome...");
        genome = GenomeGenerator.generate(1_000_000, 42, DEFAULT_MOTIFS, PLANT_COPIES);
        genomeName = "synthetic";
        MotifScanner.Job job = MotifScanner.prepare(DEFAULT_MOTIFS, true);
        long t0 = System.nanoTime();
        MotifHit[] hits = MotifScanner.run(genome, job, "aho-corasick");
        printResults(hits, job, "aho-corasick", (System.nanoTime() - t0) / 1e6, true);
        System.out.println();
        Benchmark.run(genome, DEFAULT_MOTIFS);
        System.out.println();
        Benchmark.repeatAnalysis(genome);
    }

    private static void scale(String[] args) {
        int[] sizes = { 100_000, 1_000_000, 5_000_000, 10_000_000, 20_000_000 };
        if (args.length > 1) {
            sizes = new int[args.length - 1];
            for (int i = 1; i < args.length; i++) sizes[i - 1] = Integer.parseInt(args[i]);
        }
        System.out.println("Median search time per algorithm, motifs: " + join(DEFAULT_MOTIFS));
        Benchmark.scaling(sizes, DEFAULT_MOTIFS);
    }
}
