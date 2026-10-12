# DNA Sequence Analysis Tool - Advanced String Matching in Java

## Abstract and objectives
Biological motifs (promoter boxes such as `TATAAA`, restriction sites such as `GAATTC`, start codons such as
`ATG`) must be located inside genome-scale sequences. This project searches several motifs at once, reports every
match position (overlaps included), and compares exact string-matching algorithms on the same genome.

Objectives: (1) locate motifs accurately in large DNA sequences, (2) support many motifs simultaneously,
(3) compare algorithms fairly and verify that they agree, (4) show where suffix structures help
(longest repeated region). The course rule is kept: **no `java.util.*`**; every structure is hand-built
(`IntList`, `Lcg`, `Lcp`, `Strs`, plain arrays). Only `java.io` is used, for files.

## Requirements
* JDK 11 or newer recommended. The source compiles with `--release 8`, `11` and `17` (checked); it was run on JDK 21.
* VS Code with the "Extension Pack for Java" (optional - the terminal scripts work without it).

## Run it (Windows / VS Code)
1. Unzip and open the `DNAMotifFinder` folder in VS Code (File > Open Folder - the folder that contains `src`).
2. Open `src/texthack/app/Main.java` and press **Run** above `main`, or use Run and Debug and pick
   *DNA Sequence Analysis Tool (menu)*. The menu appears in the integrated terminal.
3. Other launch entries: *Run self-tests* and *Auto demo*.

Terminal alternative (from the project folder): `run.bat` (Windows) or `./run.sh`; add `test`, `demo` or `scale`
(for example `run.bat test`). For very large genomes use more heap: `java -Xmx4g -cp out texthack.app.Main`.

## Using the menu
```
1. Load DNA Genome from FASTA        2. Generate Synthetic Genome
3. Search DNA Motifs                 4. Compare Algorithms (Benchmark)
5. Find Longest Repeated Region      6. Save Genome as FASTA        7. Exit
```
Typical demo: `2` (Enter) -> `3` (Enter, Enter) -> `4` -> `5`. Pressing Enter accepts the default shown in `[ ]`.

Search asks only for the motifs (comma separated) and the algorithm (default Aho-Corasick). Optional extras are typed
at the end of the motif line, so they never interrupt the normal flow:
`ATG,TATAAA --rc` (also search the reverse-complement strand) and `ATG,GCG --csv hits.csv` (export every position).
Reverse-complement search is off by default. Positions are **0-based start indices** (first base = 0) and are labelled as
such in every result.

### Sample (genome `TestGenome` = `ATGCGATAGCGATGCG`, file `sample/test_genome.fasta`)
```
Genome    : TestGenome (16 bases)
Algorithm : Aho-Corasick
Positions : 0-based start index (the first base of the genome is position 0)

Motif                    Count   Positions
ATG                      2       0, 11
GCG                      3       2, 8, 13
TAG                      1       6
```

## Algorithms
| Algorithm | Purpose / idea | Time |
|---|---|---|
| Naive | baseline: compare at every shift | O(nm) worst case |
| KMP | failure function (longest border) avoids re-comparing | O(n+m) |
| Rabin-Karp | double polynomial rolling hash, every hash hit is verified | O(n+m) expected |
| Z algorithm | Z-array of `P + sep + T` | O(n+m) |
| Aho-Corasick | trie + failure + dictionary links; all motifs in one pass | O(n + sum of motif lengths + matches) |
| Suffix array | prefix doubling with counting sort; query by binary search | O(n log n) build, O(m log n + occ) query |
| LCP array (Kasai) | longest common prefixes of neighbouring suffixes | O(n) |

Aho-Corasick is the natural fit for "many motifs, one genome". The suffix array is offered as a search option and is
used for repeated-region analysis; it is deliberately **not** in the benchmark table because it solves a different
problem (index once, query many times), so timing its build against a single scan would be misleading.

## Data structures
`char[]` genome; `IntList` (growable int array, merge sort) for match positions; `int[][]` goto automaton over the
alphabet A,C,G,T,N plus `int[]` failure/dictionary links; `int[]` failure table (KMP), Z-array, suffix array, LCP array.

## Bioinformatics features
Reverse-complement (minus strand) search with palindromic sites such as `GAATTC` counted once; FASTA reader/writer
(lower case folded, IUPAC codes other than ACGT become `N`, multiple records joined by one `N` so a motif cannot span two
records); synthetic genome generator with planted motifs; longest repeated region (ignores stretches of `N`).

## Benchmark method
Same genome and same prepared motifs for every algorithm. File loading and motif validation are outside the timed
region; the timed region is each algorithm's own preprocessing plus the scan (so Aho-Corasick's automaton construction is
included). A warm-up phase lets the JIT settle, then the search is repeated (up to 7 times, fewer when one run is slow) and the
**median** is reported using `System.nanoTime()`. Every algorithm's positions are compared with Naive's, and the table says so.
Times vary between runs and between computers, and no algorithm is fastest on every input.

### Measured results (my development environment, not your PC)
Linux, OpenJDK 21, `-Xmx2g`, synthetic genome with 6 motifs (forward strand). Your numbers will differ.

1,000,000 bases (19,379 matches, identical for all algorithms):

| Algorithm | Match Count | Time (ms) |
|---|---|---|
| Naive | 19379 | 31.6 |
| KMP | 19379 | 31.6 |
| Rabin-Karp | 19379 | 39.3 |
| Aho-Corasick | 19379 | 9.9 |
| Z algorithm | 19379 | 45.4 |

Increasing genome size (`scale` command, median ms, all algorithms agreed at every size):

| Bases | Naive | KMP | Rabin-Karp | Aho-Corasick | Z |
|---|---|---|---|---|---|
| 100,000 | 3.3 | 3.2 | 4.6 | 1.1 | 6.8 |
| 1,000,000 | 37.2 | 40.1 | 45.2 | 10.5 | 54.7 |
| 5,000,000 | 210.4 | 172.5 | 262.3 | 56.2 | 330.6 |
| 10,000,000 | 367.9 | 352.5 | 415.9 | 95.6 | 763.6 |
| 20,000,000 | 690.3 | 647.1 | 887.3 | 197.6 | 1131.6 |

Reading the numbers: on random DNA a mismatch usually happens within the first letters, so Naive behaves almost
linearly and is close to KMP; the guarantee of KMP/Z/Rabin-Karp/Aho-Corasick matters on repetitive input. Worst-case
demo (400,000 x `A`, motif = 599 x `A` + `C`): Naive 275.9 ms, KMP 1.3 ms, Rabin-Karp 2.6 ms, Aho-Corasick 1.4 ms, Z 4.9 ms.
Aho-Corasick leads here partly because it handles all six motifs in one pass; with a single motif the gap would be smaller.
Suffix array (1,000,000 bases): built in about 154 ms, LCP in about 37 ms, longest repeat (19 bases) found in about 13 ms.

## Testing
`Main test` (or the *Run self-tests* launch entry) runs 39 checks with no external libraries:
* the specification genome: ATG -> 0,11; GCG -> 2,8,13; TAG -> 6 (every algorithm);
* a motif that does not occur, single-character motif, overlapping matches (AAA in AAAAA -> 0,1,2), motif at the start and
  at the end, motif equal to / longer than the genome, lower-case input, multiple motifs, duplicate motifs, invalid input;
* reverse-complement search and palindromic sites; FASTA save/load round trip, header, N handling, missing file;
* 3000 random texts x 2 motifs x 2 strand modes x 6 algorithms compared with Naive (0 mismatches);
* benchmark correctness (all algorithms identical on a 200,000-base genome);
* suffix array and Kasai LCP against brute force on random strings, longest repeat (including the `N` case);
* genome sizes 1,000 / 10,000 / 100,000 / 1,000,000, all six algorithms agreeing.

Current status: **39 passed, 0 failed** (run in my environment on JDK 21).

## Implementation status
Implemented and tested: everything above. Verified only in my Linux environment: compilation, the self-tests,
scripted runs of every menu option, matching up to 20,000,000 synthetic bases, suffix-array analysis up to 1,000,000 bases.

## Limitations and future improvements
* Exact matching only; no mismatches and no IUPAC wildcard motifs (R, Y, ...). Future: mismatch-tolerant search with edit distance.
* The whole genome is held in memory as `char[]` (2 bytes per base). The menu accepts up to 50,000,000 bases, but only up to
  20,000,000 were tested; 100 MB genomes are untested and not claimed. Future: 2-bit packing, streaming.
* Suffix array uses O(n log n) prefix doubling, not SA-IS; analysing very large genomes will be slow.
* The longest-repeat result is the longest repeated substring; its two copies may overlap.
* Positions are `int` offsets into the joined multi-record sequence (no per-record coordinates).
* Not tested: Windows and `run.bat`, the VS Code launch entries on a student's machine, real genome FASTA files.

## Layout
`util/` IntList, Dna, Lcg, Lcp, Strs | `algo/` matchers, AhoCorasick, SuffixArrayIndex | `bio/` MotifScanner, MotifHit,
Fasta, GenomeGenerator | `app/` Main, Benchmark, SelfTest | `sample/` sample_genome.fasta (200 kb), test_genome.fasta (TestGenome)
