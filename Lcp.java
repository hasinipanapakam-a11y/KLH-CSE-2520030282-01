package texthack.util;

/**
 * LCP array (Kasai, O(n)) and the repeated-region search built on top of it.
 * lcp[i] = length of the longest common prefix of the suffixes sa[i-1] and sa[i] (lcp[0] = 0).
 */
public final class Lcp {
    private Lcp() {}

    /** Kasai's linear-time LCP construction from a suffix array. */
    public static int[] kasai(char[] text, int[] sa) {
        int n = text.length;
        int[] rank = new int[n];
        int[] lcp = new int[n];
        for (int i = 0; i < n; i++) rank[sa[i]] = i;
        int h = 0;
        for (int i = 0; i < n; i++) {
            if (rank[i] > 0) {
                int j = sa[rank[i] - 1];
                while (i + h < n && j + h < n && text[i + h] == text[j + h]) h++;
                lcp[rank[i]] = h;
                if (h > 0) h--;
            } else h = 0;
        }
        return lcp;
    }

    /**
     * Longest substring that occurs at least twice, ignoring any stretch containing 'N'
     * (runs of unknown bases are not a meaningful repeat). Returns {startA, startB, length};
     * length 0 means no repeat exists. The two occurrences may overlap. O(n).
     */
    public static int[] longestRepeat(char[] text, int[] sa, int[] lcp) {
        int n = text.length;
        int[] noN = new int[n + 1];                 // noN[i] = number of non-N chars starting at i
        for (int i = n - 1; i >= 0; i--) noN[i] = text[i] == 'N' ? 0 : noN[i + 1] + 1;
        int best = 0, a = 0, b = 0;
        for (int i = 1; i < n; i++) {
            int len = Math.min(lcp[i], Math.min(noN[sa[i]], noN[sa[i - 1]]));
            if (len > best) { best = len; a = sa[i - 1]; b = sa[i]; }
        }
        if (a > b) { int t = a; a = b; b = t; }
        return new int[] { a, b, best };
    }
}
