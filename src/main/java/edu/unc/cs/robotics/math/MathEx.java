package edu.unc.cs.robotics.math;

/**
 * Some possibly useful math routines
 */
public final class MathEx {
    private MathEx() {
        throw new AssertionError("static utility class, no instances");
    }

    /**
     * Casts an array of doubles to an array of floats.
     *
     * @param dst where to store results
     * @param src source of copy
     */
    public static void convert(float[] dst, double[] src) {
        final int n = dst.length;
        if (n != src.length) {
            throw new IllegalArgumentException();
        }
        for (int i=0 ; i<n ; ++i) {
            dst[i] = (float)src[i];
        }
    }

    /**
     * Casts an array of floats to an array of doubles.
     *
     * @param dst where to store results
     * @param src source of copy
     */
    public static void convert(double[] dst, float[] src) {
        final int n = dst.length;
        if (n != src.length) {
            throw new IllegalArgumentException();
        }
        for (int i=0 ; i<n ; ++i) {
            dst[i] = src[i];
        }
    }

    /**
     * Returns a copy of the array cast to floats.
     *
     * @param doubles the array of doubles to convert
     * @return an array of floats
     */
    public static float[] toFloat(double[] doubles) {
        final float[] floats = new float[doubles.length];
        convert(floats, doubles);
        return floats;
    }

    /**
     * Returns a copy of the array cast to doubles.
     *
     * @param floats the array of floats to convert
     * @return an array of floats
     */
    public static double[] toDouble(float[] floats) {
        final double[] doubles = new double[floats.length];
        convert(doubles, floats);
        return doubles;
    }

    /**
     * Tests for bit-wise equality between two floats.
     * This is computationally equivalent to:
     * {@code Float.floatToIntBits(a) == Float.floatToIntBits(b)}.
     * <p>
     * The alternate {@code a == b} fails for the cases
     * {@code 0f == -0f} (true, but should be false) and
     * {@code NaN == NaN} (false, but should be true.
     * <p>
     * The alternate {@code Float.compare(a,b) == 0} fails
     * for {@code 0f == -0f}.
     * <p>
     * The alternate fo comparing floatToIntBits directly is
     * slow.
     *
     * @param a value to compare
     * @param b value to compare
     * @return true {@code a} and {@code b} have the same bit
     * representation.
     */
    public static boolean bitEquals(float a, float b) {
        // Comparison of timings in ms for all possible floating point
        // bit representations of exactEquals(a,a):

        // return Float.floatToIntBits(a) == Float.floatToIntBits(b);
        // measured 8795.758874

        // return Float.compare(a, b) == 0;
        // measured 11586.925963 (and it fails for 0 == -0)

        // the following measured: 2197.5799859999997

        if (a == b) {
            return a != 0f || Float.floatToIntBits(a) == Float.floatToIntBits(b);
        } else {
            return Float.isNaN(a) && Float.isNaN(b);
        }
    }

    /**
     * @see #bitEquals(float, float)
     * @param a value to compare
     * @param b value to compare
     * @return true {@code a} and {@code b} have the same bit
     * representation.
     */
    public static boolean bitEquals(double a, double b) {
        // return Double.compare(a,b) == 0
        // 7025.4628729999995

        // return Double.doubleToLongBits(a) == Double.doubleToLongBits(b);
        // 6838.025492

        if (a == b) {
            return a != 0 || Double.doubleToLongBits(a) == Double.doubleToLongBits(b);
        } else {
            return Double.isNaN(a) && Double.isNaN(b);
        }
    }
}
