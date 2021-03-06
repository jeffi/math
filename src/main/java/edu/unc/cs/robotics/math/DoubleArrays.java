package edu.unc.cs.robotics.math;

import java.util.Arrays;

/**
 * Mathematical operations on arrays of double.
 */
public final class DoubleArrays {
    private DoubleArrays() {
        throw new AssertionError("no instances");
    }

    /**
     * Computes r[i] = a[i] + b[i] for all i.
     *
     * @param r result array (may be either a or b)
     * @param a left-hand side operand
     * @param b right-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static void add(double[] r, double[] a, double[] b) {
        int i = r.length;
        if (a.length != i || b.length != i) {
            throw new IllegalArgumentException();
        }
        while (--i >= 0) {
            r[i] = a[i] + b[i];
        }
    }

    /**
     * Computes r[i] = a[i] + v for all i.
     *
     * @param r result array (may be a)
     * @param a left-hand side operand
     * @param v value to add
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static void add(double[] r, double[] a, double v) {
        int i = r.length;
        if (a.length != i) {
            throw new IllegalArgumentException();
        }
        while (--i >= 0) {
            r[i] = a[i] + v;
        }
    }

    /**
     * Computes r[i] = a[i] - b[i] for all i.
     *
     * @param r result array (may be either a or b)
     * @param a left-hand side operand
     * @param b right-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static void sub(double[] r, double[] a, double[] b) {
        int i = r.length;
        if (a.length != i || b.length != i) {
            throw new IllegalArgumentException();
        }
        while (--i >= 0) {
            r[i] = a[i] - b[i];
        }
    }

    /**
     * Computes r[i] = a[i] + v for all i.
     *
     * @param r result array (may be either a or b)
     * @param a left-hand side operand
     * @param v right-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static void sub(double[] r, double[] a, double v) {
        int i = r.length;
        if (a.length != i) {
            throw new IllegalArgumentException();
        }
        while (--i >= 0) {
            r[i] = a[i] - v;
        }
    }

    /**
     * Computes r[i] = -a[i] for all i.
     *
     * @param r result array (may be either a or b)
     * @param a left-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static void neg(double[] r, double[] a) {
        int i = r.length;
        if (a.length != i) {
            throw new IllegalArgumentException();
        }
        while (--i >= 0) {
            r[i] = -a[i];
        }
    }

    /**
     * Computes the squared Euclidean distance between two vectors
     * as represented by the arrays, thus sum((b[i] - a[i])^2 for all i)
     *
     * @param a left-hand side operand
     * @param b right-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static double distSquared(double[] a, double[] b) {
        int i = a.length;
        if (b.length != i) {
            throw new IllegalArgumentException();
        }
        double sum = 0.0;
        while (--i >= 0) {
            double d = b[i] - a[i];
            sum += d*d;
        }
        return sum;
    }

    /**
     * Computes r[i] = a[i] + b[i] for all i.
     *
     * @param a left-hand side operand
     * @param b right-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static double dist(double[] a, double[] b) {
        return Math.sqrt(distSquared(a, b));
    }

    /**
     * Computes r[i] = a[i] * b[i] for all i.
     *
     * @param r result array (may be either a or b)
     * @param a left-hand side operand
     * @param b right-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static void mul(double[] r, double[] a, double[] b) {
        int i = r.length;
        if (a.length != i || b.length != i) {
            throw new IllegalArgumentException();
        }
        while (--i >= 0) {
            r[i] = a[i] * b[i];
        }
    }

    /**
     * Computes r[i] = a[i] * v for all i.
     *
     * @param r result array (may be either a or b)
     * @param a left-hand side operand
     * @param v right-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static void mul(double[] r, double[] a, double v) {
        int i = r.length;
        if (a.length != i) {
            throw new IllegalArgumentException();
        }
        while (--i >= 0) {
            r[i] = a[i] * v;
        }
    }

    /**
     * Computes r[i] = a[i] / v for all i.
     *
     * @param r result array (may be either a or b)
     * @param a left-hand side operand
     * @param v right-hand side operand
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static void div(double[] r, double[] a, double v) {
        int i = r.length;
        if (a.length != i) {
            throw new IllegalArgumentException();
        }
        while (--i >= 0) {
            r[i] = a[i] / v;
        }
    }

    /**
     * Computes r[i] = sum(a[i] * b[i] for all i).
     *
     * @param a left-hand side operand
     * @param b right-hand side operand
     * @return sum of products
     * @throws IllegalArgumentException if array sizes do not match
     */
    public static double dot(double[] a, double[] b) {
        int i = a.length;
        if (b.length != i) {
            throw new IllegalArgumentException();
        }
        double result = 0.0;
        while (--i >= 0) {
            result += a[i]*b[i];
        }
        return result;
    }

    /**
     * Computes the squared L2 length of the vector.
     *
     * @param a the vector
     * @return the squared length
     */
    public static double lengthSquared(double[] a) {
        return dot(a, a);
    }

    /**
     * Computes the length of the vector
     *
     * @param a the vector
     * @return the length of the vector
     */
    public static double length(double[] a) {
        return Math.sqrt(dot(a, a));
    }

    /**
     * Computes the maximum value contained within the argument array.
     *
     * @param values the values
     * @return the maximum value
     * @throws IndexOutOfBoundsException if values is empty
     */
    public static double max(double ... values) {
        int i = values.length;
        double r = values[--i], v;
        while (i > 0) {
            if ((v = values[--i]) > r)
                r = v;
        }
        return r;
    }

    /**
     * Computes the minimum value contained within the argument array.
     *
     * @param values the values
     * @return the minimum value
     * @throws IndexOutOfBoundsException if values is empty
     */
    public static double min(double ... values) {
        int i = values.length;
        double r = values[--i], v;
        while (i > 0) {
            if ((v = values[--i]) < r)
                r = v;
        }
        return r;
    }

    /**
     * Computes the maximum value contained within the argument array.
     *
     * @param values the values
     * @return the maximum value.  (0 if array is empty)
     */
    public static double sum(double ... values) {
        double sum = 0.0;
        for (int i=values.length ; --i >= 0 ; ) {
            sum += values[i];
        }
        return sum;
    }

    /**
     * Sets the matrix to be the identity matrix.
     *
     * @param r the matrix to set
     * @param m the number of rows and columns in r
     */
    public static void identity(double[] r, int m) {
        final int len = r.length, inc=m+1;
        if (m*m != len)
            throw new IllegalArgumentException();

        Arrays.fill(r, 0.0);
        for (int i=0 ; i<len ; i+=inc)
            r[i] = 1.0;
    }

    /**
     * Multiplies matrix a by b and stores the result in r.  This will
     * multiply any compatible pair of matrix sizes stored in column-major
     * order.
     *
     * @param r result matrix
     * @param a left-hand side operand
     * @param b right-hand side operand
     * @param m the number of rows in r and a
     * @param o the number of columns in a, and number of rows in b
     * @param n the number of columns in r and b
     */
    public static void mul(double[] r, double[] a, double[] b, int m, int o, int n) {
        if (m*n != r.length)
            throw new IllegalArgumentException();
        if (m*o != a.length)
            throw new IllegalArgumentException();
        if (o*n != b.length)
            throw new IllegalArgumentException();

        if (r != a && r != b) {
            for (int j=0, ri=0, bOff=0 ; j<n ; ++j, bOff += o) {
                for (int i=0 ; i<m ; ++i, ++ri) {
                    double s = 0.0;
                    for (int k=0, ai=i ; k<o ; ++k, ai += m)
                        s += a[ai] * b[bOff + k];
                    r[ri] = s;
                }
            }
        } else {
            // when the result matrix is one of the operands, we need
            // a temporary row of variables to avoid overwriting our
            // operands.
            final double[] bCol = new double[o];
            for (int j=0, ri=0 ; j<n ; ++j) {
                System.arraycopy(b, j*o, bCol, 0, o);
                for (int i=0 ; i<m ; ++i, ++ri) {
                    double s = 0.0;
                    for (int k=0, ai=i ; k<o ; ++k, ai += m)
                        s += a[ai] * bCol[k];
                    r[ri] = s;
                }
            }
        }
    }

    /**
     * Block matrix multiply
     *
     * @param r where to store result
     * @param rOffset offset in r
     * @param rStride stride (spacing between columns) in r
     * @param a left-hand side operand
     * @param aOffset offset in a
     * @param aStride stride (spacing between columns) in a
     * @param b right-hand side operand
     * @param bOffset offset in b
     * @param bStride stride (spacing between columns) in b
     * @param m the number of rows in r and a
     * @param o the number of columns in a, and number of rows in b
     * @param n the number of columns in r and b
     */
    public static void mul(
        double[] r, int rOffset, int rStride,
        double[] a, int aOffset, int aStride,
        double[] b, int bOffset, int bStride,
        int m, int o, int n)
    {
        if (rStride < m)
            throw new IllegalArgumentException();
        if (aStride < m)
            throw new IllegalArgumentException();
        if (bStride < o)
            throw new IllegalArgumentException();

        if (r != a && r != b) {
            for (int j=0, ri=rOffset, bOff=bOffset ; j<n ; ++j, bOff += bStride, ri += rStride - m) {
                for (int i=0 ; i<m ; ++i, ++ri) {
                    double s = 0.0;
                    for (int k=0, ai=aOffset + i ; k<o ; ++k, ai += aStride)
                        s += a[ai] * b[bOff + k];
                    r[ri] = s;
                }
            }
        } else {
            // when the result matrix is one of the operands, we need
            // a temporary row of variables to avoid overwriting our
            // operands.
            final double[] bCol = new double[o];
            for (int j=0, ri=rOffset ; j<n ; ++j, ri += rStride - m) {
                System.arraycopy(b, bOffset + j*bStride, bCol, 0, o);
                for (int i=0 ; i<m ; ++i, ++ri) {
                    double s = 0.0;
                    for (int k=0, ai=i ; k<o ; ++k, ai += aStride)
                        s += a[ai] * bCol[k];
                    r[ri] = s;
                }
            }
        }
    }
}
