package edu.unc.cs.robotics.math;

/**
 * Created by jeffi on 10/31/16.
 */
public class MathEx {
    private MathEx() {
        throw new AssertionError("static utility class, no instances");
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
