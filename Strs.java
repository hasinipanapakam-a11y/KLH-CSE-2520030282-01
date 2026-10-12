package texthack.util;

/** Small string helpers without java.util. */
public final class Strs {
    private Strs() {}

    public static String[] splitComma(String s) {
        int cnt = 1;
        for (int i = 0; i < s.length(); i++) if (s.charAt(i) == ',') cnt++;
        String[] r = new String[cnt];
        int k = 0, start = 0;
        for (int i = 0; i <= s.length(); i++) {
            if (i == s.length() || s.charAt(i) == ',') {
                r[k++] = s.substring(start, i).trim();
                start = i + 1;
            }
        }
        return r;
    }

    public static String repeat(char c, int n) {
        char[] a = new char[n];
        for (int i = 0; i < n; i++) a[i] = c;
        return new String(a);
    }
}
