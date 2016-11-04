package edu.unc.cs.robotics.math;

import java.util.Arrays;

import junit.framework.TestCase;

/**
 * Created by jeffi on 8/2/16.
 */
public class PolynomialsTest extends TestCase {
    private static void checkOne(double a, double b, double c, double d,
                              double ... expectedRoots)
    {
        double[] actualRoots = new double[3];
        int n = Polynomials.solve(actualRoots, a, b, c, d);
        assertEquals(expectedRoots.length, n);

        Arrays.sort(expectedRoots);
        Arrays.sort(actualRoots, 0, n);

        for (int i=0 ; i<n ; ++i) {
            assertEquals("root "+i, expectedRoots[i], actualRoots[i], 1e-14);
        }
    }

    private static void check(double a, double b, double c, double d,
                              double ... expectedRoots) {
        checkOne(a, b, c, d, expectedRoots.clone());

        // check that negating (e.g. reflecting about x-axis) produces
        // the same roots
        checkOne(-a, -b, -c, -d, expectedRoots.clone());

        // check that negating the odd terms (e.g. reflecting about the y-axis
        // produces the negative roots
        double[] negativeRoots = new double[expectedRoots.length];
        for (int i=0 ; i<negativeRoots.length ; ++i) {
            negativeRoots[i] = -expectedRoots[i];
        }
        checkOne(-a, b, -c, d, negativeRoots.clone());

        // check that reflecting about both axes produces negative
        // results as well.
        checkOne(a, -b, c, -d, negativeRoots.clone());
    }

    public void testSimpleCubic() throws Exception {
        check(1, 0, 0, 0, /* = */ 0.0);
    }

    public void testCubicPlusConstant() throws Exception {
        check(1, 0, 0, 27, /* = */ -3.0);
    }

    public void testSchectexExample() {
        check(1, 0, -15, -4, /* = */ 4.0, -2 - Math.sqrt(3), Math.sqrt(3) - 2);
    }

    public void testExample() throws Exception {
        check(1, 4, -8, -7, /* = */
            1.948690417436908,
            -5.266630838361576,
            -0.682059579075332);
    }

    public void testQuadratic() throws Exception {
        check(0, 2, 3, -5, /* = */ 1.0, -5.0/2.0);
    }

    public void testCubicWith2Roots() throws Exception {
        check(1, 1, 0, 0, /* = */ -1.0, 0.0);
    }

    public void testSimpleQuadratic() throws Exception {
        check(0, 1, 0, 0, /* = */ 0.0);
    }

    public void testQuadraticWithNoRoots() throws Exception {
        check(0, 1, 0, 1);
    }

    public void testLinear() throws Exception {
        check(0, 0, 3, -5, /* = */ 5.0/3.0);
    }

    public void testDegenerate() throws Exception {
        check(0, 0, 0, 0);
    }

    //    public void testBenchmark() throws Exception {
//        final int N = 1000000;
//        double sum = 0.0;
//        for (int iter = 0 ; iter<10 ; ++iter) {
//            final long t0 = System.nanoTime();
//            for (int i = 0 ; i < N ; ++i) {
//                double x = 100.0/N;
//                sum += Math.pow(x, 1.0/3.0);
//            }
//            final long t1 = System.nanoTime();
//            for (int i = 0 ; i < N ; ++i) {
//                double x = 100.0/N;
//                sum += Math.cbrt(x);
//            }
//            final long t2 = System.nanoTime();
//            System.out.printf("pow  = %.0fus\tcbrt = %.0fus (%f)%n", (t1 - t0) * 1e-3, (t2 - t1) * 1e-3, sum);
//        }
//
//    }
}