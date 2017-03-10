package edu.unc.cs.robotics.math;

import java.util.Random;

public final class MersenneTwister extends Random {

    private static final int WORD_SIZE = 64;
    private static final int STATE_SIZE = 312;
    private static final int SHIFT_SIZE = 156;
    private static final int MASK_BITS = 31;
    private static final long XOR_MASK = 0xb5026f5aa96619e9L;
    private static final int TEMPERING_U = 29;
    private static final long TEMPERING_D = 0x5555555555555555L;
    private static final int TEMPERING_S = 17;
    private static final long TEMPERING_B = 0x71d67fffeda60000L;
    private static final int TEMPERING_T = 37;
    private static final long TEMPERING_C = 0xfff7eee000000000L;
    private static final int TEMPERING_L = 43;
    private static final long INIT_MULTIPLIER = 6364136223846793005L;

    static final int word_size = WORD_SIZE;
    static final int state_size = STATE_SIZE;
    static final int shift_size = SHIFT_SIZE;
    static final int mask_bits = MASK_BITS;
    static final long xor_mask = XOR_MASK;
    static final int tempering_u = TEMPERING_U;
    static final long tempering_d = TEMPERING_D;
    static final int tempering_s = TEMPERING_S;
    static final long tempering_b = TEMPERING_B;
    static final int tempering_t = TEMPERING_T;
    static final long tempering_c = TEMPERING_C;
    static final int tempering_l = TEMPERING_L;
    static final long initialization_multiplier = INIT_MULTIPLIER;
    static final long default_seed = 5489;

    static final long parameter_a = XOR_MASK;
    static final int output_u = TEMPERING_U;
    static final int output_s = TEMPERING_S;
    static final long output_b = TEMPERING_B;
    static final int output_t = TEMPERING_T;
    static final long output_c = TEMPERING_C;
    static final int output_l = TEMPERING_L;

    static final boolean has_fixed_range = false;

    private long[] x;
    int i;

    public MersenneTwister() {
        // the super() calls setSeed which sets x.
        if (x == null)
            throw new AssertionError();
    }

    public MersenneTwister(long value) {
        super(value);

        // the super() calls setSeed which sets x.
        if (x == null)
            throw new AssertionError();
    }

    private void normalize_state() {
        final long upper_mask = -1L << MASK_BITS;
        final long lower_mask = ~upper_mask;

        long y0 = x[SHIFT_SIZE - 1] ^ x[STATE_SIZE - 1];
        if ((y0 & (1L << (WORD_SIZE - 1L))) != 0) {
            y0 = ((y0 ^ XOR_MASK) << 1) | 1;
        } else {
            y0 <<= 1;
        }

        x[0] = (x[0] & upper_mask) | (y0 & lower_mask);

        for (int j = 0 ; j < STATE_SIZE ; ++j) {
            if (x[j] != 0)
                return;
        }
        x[0] = 1L << (WORD_SIZE - 1);
    }


    private static final long UPPER_MASK = -1L << MASK_BITS;
    private static final long LOWER_MASK = ~UPPER_MASK;

    /**
     * Performs a full "twist" on the backing array.  This code is split
     * to avoid modulus for a possible performance benefit.  It is only
     * called during seeding.
     */
    private void twist() {
        for (int j = 0 ; j < STATE_SIZE - SHIFT_SIZE ; ++j) {
            long y = (x[j] & UPPER_MASK) | (x[j + 1] & LOWER_MASK);
            x[j] = x[j + SHIFT_SIZE] ^ (y >>> 1) ^ ((x[j + 1] & 1)*XOR_MASK);
        }
        for (int j = STATE_SIZE - SHIFT_SIZE ; j < STATE_SIZE - 1; ++j) {
            long y = (x[j] & UPPER_MASK) | (x[j + 1] & LOWER_MASK);
            x[j] = x[j - (STATE_SIZE - SHIFT_SIZE)] ^ (y >>> 1) ^ ((x[j + 1] & 1)*XOR_MASK);
        }
        // last iteration
        long y = (x[STATE_SIZE - 1] & UPPER_MASK) | (x[0] & LOWER_MASK);
        x[STATE_SIZE - 1] = x[SHIFT_SIZE - 1] ^ (y >>> 1) ^ ((x[0] & 1)*XOR_MASK);
        i = 0;
    }

    @Override
    public final long nextLong() {
        // The alternate way to do this is to twist() when out of bits
        // if (i == n) twist();

        // Our method twists one element at a time as we go.  This avoids
        // sporadic performance of the random number generation, which
        // could be problematic when trying to achieve a timing objective.
        // The downside is that it means we have to perform two modulus
        // operations per call.
        long v = x[i];
        final int j = (i + 1)%STATE_SIZE;
        final long y = (v & UPPER_MASK) | (x[j] & LOWER_MASK);
        x[i] = x[(i + SHIFT_SIZE)%STATE_SIZE] ^ (y >>> 1) ^ ((x[j] & 1)*XOR_MASK);

        i = j;

        v ^= ((v >>> TEMPERING_U) & TEMPERING_D);
        v ^= ((v << TEMPERING_S) & TEMPERING_B);
        v ^= ((v << TEMPERING_T) & TEMPERING_C);
        v ^= (v >>> TEMPERING_L);

        return v;
    }

    @Override
    public void setSeed(long seed) {
        // initialization issue...
        // setSeed is called by java.util.Random() constructor
        // which is called before x is initialized
        if (x == null)
            x = new long[STATE_SIZE];

        x[0] = seed;
        for (i = 1; i < STATE_SIZE ; ++i) {
            x[i] = (INIT_MULTIPLIER*(x[i - 1] ^ (x[i - 1] >>> (WORD_SIZE - 2))) + i);
        }

        normalize_state();
        twist();
    }

    @Override
    public void nextBytes(byte[] bytes) {
        int j;
        for (j = 0; j + 8 < bytes.length ; j += 8) {
            long v = nextLong();
            bytes[j] = (byte)(v >>> 56);
            bytes[j + 1] = (byte)(v >>> 48);
            bytes[j + 2] = (byte)(v >>> 40);
            bytes[j + 3] = (byte)(v >>> 32);
            bytes[j + 4] = (byte)(v >>> 24);
            bytes[j + 5] = (byte)(v >>> 16);
            bytes[j + 6] = (byte)(v >>> 8);
            bytes[j + 7] = (byte)v;
        }
        if (j < bytes.length) {
            long v = nextLong();
            switch (bytes.length - j) {
            default:
                throw new AssertionError();
            case 7:
                bytes[j + 6] = (byte)(v >>> 8);
            case 6:
                bytes[j + 5] = (byte)(v >>> 16);
            case 5:
                bytes[j + 4] = (byte)(v >>> 24);
            case 4:
                bytes[j + 3] = (byte)(v >>> 32);
            case 3:
                bytes[j + 2] = (byte)(v >>> 40);
            case 2:
                bytes[j + 1] = (byte)(v >>> 48);
            case 1:
                bytes[j] = (byte)(v >>> 56);
            }
        }
    }

    @Override
    protected int next(int bits) {
        return (int)(nextLong() & ~(-1 << bits));
    }

    @Override
    public int nextInt() {
        return (int)nextLong();
    }

    @Override
    public double nextDouble() {
        // the Java way

        return (nextLong() & ~(-1L << 53)) * 0x1.0p-53;

        // Boost's way (approximated)
//        long v = nextLong();
//        return v < 0
//            ? (v >>> 1) * 0x1.0p-63
//            : v * 0x1.0p-64;
    }
}
