package texthack.util;

/** Hand-built growable int array (no java.util). Used to store match positions. */
public final class IntList {
    private int[] a = new int[16];
    private int n;

    public void add(int v) {
        if (n == a.length) {
            int[] b = new int[a.length * 2];
            for (int i = 0; i < n; i++) b[i] = a[i];
            a = b;
        }
        a[n++] = v;
    }

    public int get(int i) { return a[i]; }
    public int size() { return n; }

    /** Merge sort, O(n log n). */
    public void sort() {
        if (n < 2) return;
        mergeSort(a, new int[n], 0, n);
    }

    private static void mergeSort(int[] x, int[] b, int lo, int hi) {
        if (hi - lo < 2) return;
        int mid = (lo + hi) >>> 1;
        mergeSort(x, b, lo, mid);
        mergeSort(x, b, mid, hi);
        int i = lo, j = mid, k = lo;
        while (i < mid && j < hi) b[k++] = x[i] <= x[j] ? x[i++] : x[j++];
        while (i < mid) b[k++] = x[i++];
        while (j < hi) b[k++] = x[j++];
        for (k = lo; k < hi; k++) x[k] = b[k];
    }

    public boolean sameAs(IntList o) {
        if (o.n != n) return false;
        for (int i = 0; i < n; i++) if (a[i] != o.a[i]) return false;
        return true;
    }
}
