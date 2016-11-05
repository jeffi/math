// Automatically generated from DoubleArraysTest.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/5/16.
 */
public class FloatArraysTest extends TestCase {
    private static void checkMul(
        float[] a, float[] b, int m, int o, int n,
        float ... expected)
    {
        float[] r = new float[expected.length];
        FloatArrays.mul(r, a, b, m, o, n);
        assertArrays(r, expected);
        FloatArrays.mul(r, 0, m, a, 0, m, b, 0, o, m, o, n);
        assertArrays(r, expected);

        if (a.length == r.length) {
            // test in-place overwrite
            r = a.clone();
            FloatArrays.mul(r, r, b, m, o, n);
            assertArrays(r, expected);
            r = a.clone();
            FloatArrays.mul(r, 0, m, r, 0, m, b, 0, o, m, o, n);
            assertArrays(r, expected);
        }
        if (b.length == r.length) {
            // test in-place overwrite
            r = b.clone();
            FloatArrays.mul(r, a, r, m, o, n);
            assertArrays(r, expected);
            r = b.clone();
            FloatArrays.mul(r, 0, m, a, 0, m, r, 0, o, m, o, n);
            assertArrays(r, expected);
        }
    }

    private static void assertArrays(float[] actual, float[] expected) {
        for (int i = 0 ; i < expected.length ; ++i) {
            assertEquals(expected[i], actual[i], 0.0f);
        }
    }

    public void testMatrixMul2x3x2() throws Exception {
        // | 1, 2, 3 |   |  7  8 |   |  58  64 |
        // | 4, 5, 6 | x |  9 10 | = | 139 154 |
        //               | 11 12 |
        float[] a = {
            1, 4,
            2, 5,
            3, 6
        };
        float[] b = {
            7, 9, 11,
            8, 10, 12,
        };

        checkMul(a, b, 2, 3, 2,
            58.0f,
            139.0f,
            64.0f,
            154.0f);
    }

    public void testMatrix2x2() throws Exception {
        float[] a = {
            1, 3,
            2, 4
        };
        float[] b = {
            2, 1,
            0, 2
        };

        checkMul(a, b, 2, 2, 2,
            4.0f,
            10.0f,
            4.0f,
            8.0f);
    }

    public void testMatrixMul1x3x4() throws Exception {
        //             | 13 9 7 15 |
        // [ 3 4 2 ] x |  8 7 4  6 | = [ 83 63 37 75 ]
        //             |  6 4 0  3 |
        float[] a = { 3, 4, 2 };
        float[] b = {
            13, 8, 6,
            9, 7, 4,
            7, 4, 0,
            15, 6, 3
        };
        checkMul(a, b, 1, 3, 4,
            83.0f,
            63.0f,
            37.0f,
            75.0f);
    }
}