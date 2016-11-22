package edu.unc.cs.robotics.math;

import static java.lang.Math.cbrt;
import static java.lang.Math.sqrt;

/**
 * Polynomial solver
 */
public class Polynomials {
    private static final double EPSILON = 1e-15; // Math.ulp(1.0);

    private Polynomials() {}

    private static boolean isZero(double x) {
        return Math.abs(x) <= EPSILON;
    }

    private static Complex2d R(double v) {
        return new Complex2d(v);
    }

    /**
     * Computes the roots of a quartic polynomial: {@code ax^4 + bx^3 + cx^2 + dx + e}.
     *
     * @param roots where the roots are stored upon return.  Thus array should have
     * at least 4 elements.
     * @param a the quartic term
     * @param b the cubic term
     * @param c the quadratic term
     * @param d the linear term
     * @param e the constant term
     * @return the number of roots stored in {@code roots}, possibly 0, up to 4.
     */
    public static int solve(double[] roots, double a, double b, double c, double d, double e) {
        if (isZero(a)) {
            // quartic term is missing, solve as a cubic
            return solve(roots, b, c, d, e);
        }

        // normalize to solve x^4 + bx^3 + cx^2 + dx + e
//        b /= a;
//        c /= a;
//        d /= a;
//        e /= a;

//        final double p = (8.0 * a*c - 3.0 * b*b) / (8.0 *a*a);
//        final double q = (b*b*b - 4.0*a*b*c + 8.0*a*a*d) / (8.0*a*a*a);
        final double d1 = 2.0*c*c*c - 9.0*b*c*d + 27.0*a*d*d + 27.0*b*b*e - 72.0*a*c*e;
        final double d0 = c*c - 3.0*b*d + 12.0*a*e;
        final double desc = d1*d1 - 4.0*d0*d0*d0; // discriminant * -27
        final double P = 8*a*c - 3*b*b;

        final double D = 64.0*a*a*a*e - 16.0*a*a*c*c + 16.0*a*b*b*c - 16.0*a*a*b*d - 3*b*b*b*b;
        final double Q = b*b*b - 4*a*b*c + 8*d*a*a;

//        System.out.println(desc);
        if (isZero(desc)) {
//            if (P < 0 && D < 0 && !isZero(d0)) {
//                System.out.println("CASE 3a: real double, 2 real simple");
//            } else if (D > 0 || (P > 0 && (D != 0 || Q != 0))) {
//                System.out.println("CASE 3b: real double, 2 complex");
//            } else if (d0 == 0 && D != 0) {
//                System.out.println("CASE 3c: triple root, 1 simple");
//            } else
            if (isZero(D)) {
                if (isZero(d0)) {
//                    System.out.println("CASE 3f: 4 roots -b/4a");
                    roots[0] = -b/(4.0*a);
                    return 1;
//                } else if (P < 0) {
//                    System.out.println("CASE 3d: 2 real double roots");
//                } else if (P > 0 && Q == 0.0) {
//                    System.out.println("CASE 3e: 2 complex double roots");
//                } else {
//                    System.out.println("CASE 3?: ???");
                }
            }
//        } else if (desc < 0) {
//            System.out.println("CASE 1: 2 real, 2 complex");
//        } else { // desc > 0
//            if (P < 0 || D < 0) {
//                System.out.println("CASE 2a: 4 distinct real roots");
//            } else {
//                System.out.println("CASE 2b: 2 pairs of non-real");
//            }
        }

        Complex2d[] Qs = new Complex2d(desc).sqrt().add(d1).div(2.0).cbrt();

        Complex2d q = Qs[2];
        Complex2d p3 = R(d0).div(q.clone().mul(3*a) ).add(q.clone().div(3.0*a));
        double v = P/(-12.0*a*a);
        double w = (Q)/(-4*a*a*a);

        Complex2d p4 = p3.clone().add(v).sqrt();
        Complex2d p5 = R(v*2).sub(p3);
        Complex2d p6 = R(w).div(p4);
        p4.div(2);


        Complex2d r0 = R(b/(-4*a)).sub(p4).sub(p5.clone().sub(p6).sqrt().div(2));
        Complex2d r1 = R(b/(-4*a)).sub(p4).add(p5.clone().sub(p6).sqrt().div(2));
        Complex2d r2 = R(b/(-4*a)).add(p4).sub(p5.clone().add(p6).sqrt().div(2));
        Complex2d r3 = R(b/(-4*a)).add(p4).add(p5.clone().add(p6).sqrt().div(2));


            int n = 0;
            if (isZero(r0.imaginary)) {
                roots[n++] = r0.real;
            }
            if (isZero(r1.imaginary)) {
                roots[n++] = r1.real;
            }
            if (isZero(r2.imaginary)) {
                roots[n++] = r2.real;
            }
            if (isZero(r3.imaginary)) {
                roots[n++] = r3.real;
            }


        // remove duplicates
        for (int i=0 ; i<n ; ++i) {
            for (int j=i+1 ; j<n ; ) {
                if (Math.abs(roots[i] - roots[j]) < EPSILON*4) {
                    roots[j] = roots[--n];
//                    System.out.println("removed duplicate root: "+roots[i]);
                } else {
                    ++j;
                }
            }
        }

        // insertion sort
        for (int i=0, j ; (j=i)<n-1 ; ) {
            final double ai = roots[++i];
            while (ai < roots[j]) {
                roots[j+1] = roots[j];
                if (j-- == 0) {
                    break;
                }
            }
            roots[j+1] = ai;
        }

        return n;
    }

    /**
     * Computes the roots of ax^3 + bx^2 + cx + d = 0.
     *
     * @param roots where the roots are stored upon return.  This
     * array should have at least 3 elements.
     * @param a cubic term
     * @param b quadratic term
     * @param c linear term
     * @param d constant term
     * @return number of roots stored in {@code roots}, possibly 0, up to 3.
     */
    public static int solve(double[] roots, double a, double b, double c, double d) {
        // Cubic equation from:
        // http://www.math.vanderbilt.edu/~schectex/courses/cubic/

        if (isZero(a)) {
            // cubic term is missing, solve using quadratic formula
            return solve(roots, b, c, d);
        }

        final double p = -b/(3.0*a);
        final double q = p*p*p + (b*c - 3.0*a*d)/(6.0*a*a);
        final double r = c/(3.0*a);

        final double s = r - p*p;
        final double t = q*q + s*s*s;

        if (isZero(t)) {
            // When t == 0, there are a maximum of 2 roots.  One that pass through the axis
            // and another that just grazes it with a derivative of 0.

            // This is a simplification of the complex math case
            // in which atan(im,q) == 0, thus we get angles 0, 2pi/3, 4pi/3
            // cos(0) = 1, sin(0) = 0, so we just use the cbrt
            // cos(2pi/3) = -0.5, sin(2pi/3) =  0.866...
            // cos(4pi/3) = -0.5, sin(4pi/3) = -0.866...
            // Thus for angle 0, we just add cbrt+cbrt+p
            // For the other root, we add the result of 2pi/3 and 4pi/3 since
            // the imaginary parts cancel.  If we were to run through
            // the imaginary banch code, we'd get that answer twice.

            if (isZero(q)) {
                // in this case the root that grazes is also the root
                // that passes through.
                roots[0] = p;
                return 1;
            } else {
                final double cbrt = cbrt(q);

                // Return roots in increasing order
                if (q < 0) {
                    roots[0] = 2.0*cbrt + p;
                    roots[1] = p - cbrt;
                } else {
                    roots[0] = p - cbrt;
                    roots[1] = 2.0*cbrt + p;
                }

                // (1 + 0 i)*2 + p
                // (-0.5 + 0.866 i) + (-0.5 - 0.866 i) + p
                return 2;
            }
        } else if (t > 0) {
            // if t >= 0 we can use real math (not complex),
            // and there is only one root
            double sqrt = sqrt(t);
            double c1 = cbrt(q + sqrt);
            double c2 = cbrt(q - sqrt);
            roots[0] = c1 + c2 + p;
            return 1;
        } else {
            // if t < 0, then we need to use complex math
            // Since we need sqrt(t)
            //  The following code implements
            // cbrt(q + sqrt(t)) + cbrt(q - sqrt(t)) + p
            // Since cbrt on a complex number has up to 3 complex roots
            // we only consider the roots that when summed have real
            // results.
            // Using Apache Commons for Complex, the code looks like this:

//                System.out.println("t=" + t);
//                Complex sqrt = new Complex(t, 0).sqrt();
//                System.out.println("sqrt(t)=" + sqrt);
//
//                List<Complex> roots1 = sqrt.add(q).nthRoot(3);
//                List<Complex> roots2 = new Complex(q).subtract(sqrt).nthRoot(3);
//                System.out.println(roots1);
//                System.out.println(roots2);
//
//                System.out.println("Roots====");
//                int i = 0;
//                for (Complex r1 : roots1) {
//                    for (Complex r2 : roots2) {
//                        Complex sum = r1.add(r2).add(p);
//                        if (Math.abs(sum.getImaginary()) < 1e-14) {
//                            System.out.println(sum.getReal());
//                            roots[i++] = sum.getReal();
//                        }
//                    }
//                }
//
            // cbrt(q + im)

            //==========
//                final double im = Math.sqrt(-t);
//                final double cbrtOfAbs = Math.pow(complexAbs(q, im), 1.0/3.0);
//                final double nthPhi = Math.atan2(im, q)/3.0;
//                final double slice = 2.0 * Math.PI / 3.0;
//
//                roots[0] = 2.0 * cbrtOfAbs * Math.cos(nthPhi) + p;
//                roots[1] = 2.0 * cbrtOfAbs * Math.cos(nthPhi + slice) + p;
//                roots[2] = 2.0 * cbrtOfAbs * Math.cos(nthPhi + slice*2) + p;

            final double im = sqrt(-t);
//                final double cbrtOfAbs = Math.pow(complexAbs(q, im), 1.0/3.0);
            final double cbrt = 2.0*Math.pow(q*q - t, 1.0/6.0);
            final double phi = Math.atan2(im, q)/3.0;

            // Math.cos(nthPhi) = im / Math.sqrt(im*im + q*q);


            // This ordering roots in increasing order.
            // `cbrt` is always positive since q*q is positive and (-t) is positive,
            // thus we do not have to worry about the order reversing.
            // Since `im` is also positive, atan2(im,q)/3 results in an angle that means:
            //    cos(phi+2*pi/3) < cos(phi + 4*pi/3) < cos(phi)
            // To see a graphical depiction of this, use gnuplot:
            // set yrange [0:10]
            // splot cos(atan2(y,x)/3.0 + 2*pi/3), cos(atan2(y,x)/3 - 2.0*pi/3.0), cos(atan2(y,x)/3.0)
            roots[0] = cbrt*Math.cos(phi + 2.0*Math.PI/3.0) + p;
            roots[1] = cbrt*Math.cos(phi - 2.0*Math.PI/3.0) + p;
            roots[2] = cbrt*Math.cos(phi) + p;


            // TODO: remove this assert (or move it to the tests)
            assert roots[0] < roots[1] && roots[1] < roots[2]
                : "root order incorrect: " + roots[0] + ", " + roots[1] + ", " + roots[2] +
                  " (im=" + im + ", q=" + q + ", phi=" + phi + ", cbrt="+cbrt+')';

            return 3;
        }
    }

    /**
     * Solves a for the roots of quadratic equation ax^2 + bx + c = 0.
     *
     * @param roots where the roots will be store (must be at least length 2)
     * @param a the quadratic term
     * @param b the linear term
     * @param c the constant term
     * @return the number of roots
     */
    public static int solve(double[] roots, double a, double b, double c) {
        if (isZero(a)) {
            // quadratic terms is missing, solve the
            // simple line equation.
            return solve(roots, b, c);
        }

        double term = b*b - 4.0*a*c;
        double inv2a = -0.5/a;
        if (term < 0) {
            // no solution.  parabola either completely above
            // or completely below axis
            return 0;
        } else if (term == 0.0) {
            // one solution.
            roots[0] = b*inv2a;
            return 1;
        }
        double root = sqrt(term);
        roots[0] = (b - root)*inv2a;
        roots[1] = (b + root)*inv2a;
        return 2;
    }

    /**
     * Solves for the roots of a linear equation ax + b = 0
     *
     * @param roots where the roots will be stored (must be at least length 1)
     * @param a the linear term
     * @param b the constant term
     * @return the number of roots
     */
    public static int solve(double[] roots, double a, double b) {
        if (isZero(a)) {
            return 0;
        }

        roots[0] = - b/a;
        return 1;
    }

}
