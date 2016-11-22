package edu.unc.cs.robotics.math;

import java.util.Arrays;
import java.util.Random;

import junit.framework.TestCase;

/**
 * Created by jeffi on 8/2/16.
 */
public class PolynomialsTest extends TestCase {
    private static double[] eq(double ... values) {
        return values;
    }

    private static void checkOne(
        double a, double b, double c, double d,
        double[] expectedRoots)
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

    private static void check(
        double a, double b, double c, double d,
        double[] expectedRoots)
    {
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
        check(1, 0, 0, 0, eq(0.0));
    }

    public void testCubicPlusConstant() throws Exception {
        check(1, 0, 0, 27, eq(-3.0));
    }

    public void testSchectexExample() {
        check(1, 0, -15, -4, eq(4.0, -2 - Math.sqrt(3), Math.sqrt(3) - 2));
    }

    public void testExample() throws Exception {
        check(1, 4, -8, -7, eq(
            1.948690417436908,
            -5.266630838361576,
            -0.682059579075332));
    }

    public void testQuadratic() throws Exception {
        check(0, 2, 3, -5, eq(1.0, -5.0/2.0));
    }

    public void testCubicWith2Roots() throws Exception {
        check(1, 1, 0, 0, eq(-1.0, 0.0));
    }

    public void testSimpleQuadratic() throws Exception {
        check(0, 1, 0, 0, eq(0.0));
    }

    public void testQuadraticWithNoRoots() throws Exception {
        check(0, 1, 0, 1, eq());
    }

    public void testLinear() throws Exception {
        check(0, 0, 3, -5, eq(5.0/3.0));
    }

    public void testDegenerate() throws Exception {
        check(0, 0, 0, 0, eq());
    }

    static void check(double a, double b, double c, double d, double e, double[] expected) {
        final double[] roots = new double[4];
        final int n = Polynomials.solve(roots, a, b, c, d, e);
        // System.out.println("Roots = "+Arrays.toString(Arrays.copyOf(roots, n)));
        for (int i=0 ; i<n ; ++i) {
            final double x = roots[i];
            double v = (((a*x + b)*x + c)*x + d)*x + e;
            assertEquals("root " + i + " = "+roots[i], 0.0, v, 1e-9);
        }
        assertEquals("number of roots: "+Arrays.toString(Arrays.copyOf(roots, n)), expected.length, n);
        for (int i=0 ; i<n ; ++i) {
            assertEquals(expected[i], roots[i], 1e-9);
        }
    }

    public void testSimpleQuartic() throws Exception {
        final double[] roots = new double[4];
        final int n = Polynomials.solve(roots, 1.0, 0.0, 0.0, 0.0, 0.0);
        assertEquals(1, n);
        assertEquals(0.0, roots[0], 0.0);
    }

    public void testPotatoJacksQuartic() throws Exception {
        // http://planetmath.org/comment/4497#comment-4497
        check(1.0, -2.0, 4.0, -3.0, -10.0, eq(-1.0, 2.0));
        // also has complex roots at 1/2 (1 +/- i sqrt(19))
    }

    public void testQuartic003() throws Exception {
        // x^4 - 4*x^2 = x^2 * (x^2 - 4)
        check(1.0, 0.0, -4.0, 0.0, 0.0, eq(-2.0, 0.0, 2.0));
    }

    public void testType1() throws Exception {
        check(1, -4, -2, 12, -3,
            eq(
                -1.732050807568877, // -Math.sqrt(3),
                0.26794919243112303, // 2 - Math.sqrt(3),
                1.732050807568877,  // Math.sqrt(3),
                3.7320508075688767)); // 2 + Math.sqrt(3));
    }

    public void testType2() throws Exception {
        check(1, 2, -12, -2, 6, eq(
            -4.50136748673309,
            -0.762473005315555,
            +0.67536437029438,
            +2.58847612175221));
    }

    public void testType3() throws Exception {
        // 3*x^4 - 20*x^3 + 42*x^2 - 36*x -20
        check(3, -20, 42, -36, -20,
            eq(
                -0.368164612948320,
                4.03594806101948));
    }

    public void testType4() throws Exception {
        // x^4 -6*x^2 -12*x +5
        check(1, 0, -6, -12, 5,
            eq(
                0.354983353019874,
                3.06339980193731));
    }

    public void testType5() throws Exception {
        // x^4 + 4*x^3 +6*x^2 -28*x -40
        check(1, 4, 6, -28, -40,
            eq(
                -1.28105500865338,
                2.26350273648182));
//            -2.35303360744636,
//            1.90819686194685);
    }

    public void testType6() throws Exception {
        double[] roots = new double[4];
        // x^4 -4*x^3 + 6*x^2 -4*x -3
        check(1, -4, 6, -4, -3,
            eq(
                1.0 - Math.sqrt(2),
                1.0 + Math.sqrt(2)));
    }

    public void testType7() throws Exception {
        double[] roots = new double[4];
        // x^4 - 2*x^3 + 6*x^2 + 22*x -10
        check(1, -2, 6, 22, -10,
            eq(
                -1.96461832990376,
                0.413090872794899));
    }

    public void testRandomQuartic() throws Exception {
        final Random rng = new Random(1);
        for (int i=0 ; i<100000 ; ++i) {
            double a = (rng.nextDouble() - 0.5) * 100;
            if (Math.abs(a) < 2.0)
                continue;
            double b = (rng.nextDouble() - 0.5) * 100;
            double c = (rng.nextDouble() - 0.5) * 100;
            double d = (rng.nextDouble() - 0.5) * 100;
            double e = (rng.nextDouble() - 0.5) * 100;
            double[] roots = new double[4];
            final int n = Polynomials.solve(roots, a, b, c, d, e);
            for (int j=0 ; j<n ; ++j) {
                final double x = roots[j];
                final double x2 = x*x;
                assertEquals(
                    "random #"+i+": "+a+","+b+","+c+","+d+","+e,
                    // a+"x^4 + "+b+"x^3 + "+c+"x^2 + "+d+"x + "+e,
                    0.0, a*x2*x2 + b*x2*x + c*x2 + d*x + e, 1e-9);
            }
        }
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