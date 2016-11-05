package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/5/16.
 */
public class DoubleArraysTest extends TestCase {
    private static void checkMul(
        double[] a, double[] b, int m, int o, int n,
        double ... expected)
    {
        double[] r = new double[expected.length];
        DoubleArrays.mul(r, a, b, m, o, n);
        assertArrays(r, expected);
        DoubleArrays.mul(r, 0, m, a, 0, m, b, 0, o, m, o, n);
        assertArrays(r, expected);

        if (a.length == r.length) {
            // test in-place overwrite
            r = a.clone();
            DoubleArrays.mul(r, r, b, m, o, n);
            assertArrays(r, expected);
            r = a.clone();
            DoubleArrays.mul(r, 0, m, r, 0, m, b, 0, o, m, o, n);
            assertArrays(r, expected);
        }
        if (b.length == r.length) {
            // test in-place overwrite
            r = b.clone();
            DoubleArrays.mul(r, a, r, m, o, n);
            assertArrays(r, expected);
            r = b.clone();
            DoubleArrays.mul(r, 0, m, a, 0, m, r, 0, o, m, o, n);
            assertArrays(r, expected);
        }
    }

    private static void assertArrays(double[] actual, double[] expected) {
        for (int i = 0 ; i < expected.length ; ++i) {
            assertEquals(expected[i], actual[i], 0.0);
        }
    }

    public void testMatrixMul2x3x2() throws Exception {
        // | 1, 2, 3 |   |  7  8 |   |  58  64 |
        // | 4, 5, 6 | x |  9 10 | = | 139 154 |
        //               | 11 12 |
        double[] a = {
            1, 4,
            2, 5,
            3, 6
        };
        double[] b = {
            7, 9, 11,
            8, 10, 12,
        };

        checkMul(a, b, 2, 3, 2,
            58.0,
            139.0,
            64.0,
            154.0);
    }

    public void testMatrix2x2() throws Exception {
        double[] a = {
            1, 3,
            2, 4
        };
        double[] b = {
            2, 1,
            0, 2
        };

        checkMul(a, b, 2, 2, 2,
            4.0,
            10.0,
            4.0,
            8.0);
    }

    public void testMatrixMul1x3x4() throws Exception {
        //             | 13 9 7 15 |
        // [ 3 4 2 ] x |  8 7 4  6 | = [ 83 63 37 75 ]
        //             |  6 4 0  3 |
        double[] a = { 3, 4, 2 };
        double[] b = {
            13, 8, 6,
            9, 7, 4,
            7, 4, 0,
            15, 6, 3
        };
        checkMul(a, b, 1, 3, 4,
            83.0,
            63.0,
            37.0,
            75.0);
    }
}