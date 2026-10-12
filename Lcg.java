package texthack.util;

/** Tiny deterministic random generator (SplitMix64) - replaces java.util.Random. */
public final class Lcg {
    private long state;

    public Lcg(long seed) { state = seed; }

    public long nextLong() {
        long z = (state += 0x9E3779B97F4A7C15L);
        z = (z ^ (z >>> 30)) * 0xBF58476D1CE4E5B9L;
        z = (z ^ (z >>> 27)) * 0x94D049BB133111EBL;
        return z ^ (z >>> 31);
    }

    public int nextInt(int bound) { return (int) ((nextLong() >>> 1) % bound); }
}
