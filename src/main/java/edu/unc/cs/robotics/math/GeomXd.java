package edu.unc.cs.robotics.math;

/**
 * Geometry methods that work for arbitrary dimensions (or at least 2 and 3),
 * and operate on double[].  This should be deprecated.
 */
public final class GeomXd {

    public static final double EPSILON = 1e-9;

    private GeomXd() {
        throw new AssertionError("no instances");
    }

    /**
     * Compute the time of CPA for two "tracks"
     *
     * <p>source: http://geomalgorithms.com/a07-_distance.html</p>
     *
     * <p>Track 1 starts at p1 and has velocity v1</p>
     * <p>Track 2 starts at p2 and has velocity v2</p>
     */
    public static double cpaTime(double[] p1, double[] v1,
                                 double[] p2, double[] v2)
    {
        // t_cpa = -w_0 * (u - v) / |u - v|^2
        // w_0 = P_0 - Q_0

        double dv2 = 0;
        double dot = 0;

        for (int i = p1.length ; --i >= 0 ; ) {
            final double dvi = v1[i] - v2[i];
            dv2 += dvi*dvi;
            final double w0i = p2[i] - p1[i];
            dot += w0i*dvi;
        }

        return dv2 < EPSILON ? 0 : dot/dv2;
    }

    static class DistSegmentSegment {
        final int dim;
        double[] u;
        double[] v;
        double[] w;

        DistSegmentSegment(int dim) {
            this.dim = dim;
            u = new double[dim];
            v = new double[dim];
            w = new double[dim];
        }

        double compute(double[] s1p0, double[] s1p1,
                       double[] s2p0, double[] s2p1)
        {
            DoubleArrays.sub(u, s1p1, s1p0);
            DoubleArrays.sub(v, s2p1, s2p0);
            DoubleArrays.sub(w, s1p0, s2p0);

            double a = DoubleArrays.dot(u, u);
            double b = DoubleArrays.dot(u, v);
            double c = DoubleArrays.dot(v, v);
            double d = DoubleArrays.dot(u, w);
            double e = DoubleArrays.dot(v, w);
            double D = a*c - b*b;
            double sD = D;
            double tD = D;

            double sN, tN;

            if (D < 1e-9) {
                sN = 0.0;
                sD = 1.0;
                tN = e;
                tD = c;
            } else {
                sN = (b*e - c*d);
                tN = (a*e - b*d);
                if (sN < 0) {
                    sN = 0;
                    tN = e;
                    tD = c;
                } else if (sN > sD) {
                    sN = sD;
                    tN = e + b;
                    tD = c;
                }
            }

            if (tN < 0) {
                tN = 0;
                if (-d < 0) {
                    sN = 0;
                } else if (-d > a) {
                    sN = sD;
                } else {
                    sN = -d;
                    sD = a;
                }
            } else if (tN > tD) {
                tN = tD;
                if ((-d + b) < 0.0) {
                    sN = 0;
                } else if ((-d + b) > a) {
                    sN = sD;
                } else {
                    sN = (-d + b);
                    sD = a;
                }
            }

            double sc = (Math.abs(sN) < 1e-9 ? 0.0 : sN / sD);
            double tc = (Math.abs(tN) < 1e-9 ? 0.0 : tN / tD);

            double sum = 0;
            for (int i=0 ; i<dim ; ++i) {
                double dPi = w[i] + sc * u[i] - tc * v[i];
                sum += dPi*dPi;
            }

            return Math.sqrt(sum);
        }
    }

    public static double distPointSegment(double[] pt, double[] s0, double[] s1) {
        final int dim = pt.length;

        double[] v = new double[dim];
        DoubleArrays.sub(v, s1, s0);
        double[] w = new double[dim];
        DoubleArrays.sub(w, pt, s0);
        double c1 = DoubleArrays.dot(w, v);
        double c2;
        if (c1 <= 0) {
            return DoubleArrays.length(w);
        } else if ((c2 = DoubleArrays.dot(v, v)) <= c1) {
            return DoubleArrays.dist(pt, s1);
        } else {
            double b = c1 / c2;
            double sum = 0;
            for (int i = dim; --i >= 0; ) {
                double di = s0[i] + v[i] * b - pt[i];
                sum += di * di;
            }
            return Math.sqrt(sum);
        }
    }

    /**
     * Computes the distance between a point and a line segment.
     *
     * @param pt the point
     * @param s0 endpoint of line segment
     * @param s1 endpoint of line segment
     * @param nearest on return contains the point on the segment nearest to pt.
     * @return the distance between nearest and pt.
     */
    public static double distPointSegmentSquared(double[] pt, double[] s0, double[] s1, double[] nearest) {
        final int dim = pt.length;
        double c1 = 0;
        double c2 = 0;
        // double ww = 0;
        for (int i=dim ; --i >= 0 ; ) {
            double vi = s1[i] - s0[i];
            double wi = pt[i] - s0[i];
            c1 += vi * wi;
            c2 += vi * vi;
            // ww += wi * wi;
            // borrow nearest for temporary storage.
            nearest[i] = vi;
        }
        if (c1 <= 0) {
            System.arraycopy(s0, 0, nearest, 0, dim);
            return DoubleArrays.distSquared(pt, s0);
        }
        if (c2 <= c1) {
            System.arraycopy(s1, 0, nearest, 0, dim);
            return DoubleArrays.distSquared(pt, s1);
        }
        double b = c1 / c2;
        double sum = 0;
        for (int i=dim ; --i >= 0 ; ) {
            nearest[i] = s0[i] + nearest[i] * b;
            double di = nearest[i] - pt[i];
            sum += di*di;
        }
        return sum;
    }

    public static double distSegmentSegment(
        double[] s1p0, double[] s1p1,
        double[] s2p0, double[] s2p1)
    {
        final int dim = s1p0.length;
        double[] u = new double[dim];
        double[] v = new double[dim];
        double[] w = new double[dim];
        DoubleArrays.sub(u, s1p1, s1p0);
        DoubleArrays.sub(v, s2p1, s2p0);
        DoubleArrays.sub(w, s1p0, s2p0);

        double a = DoubleArrays.dot(u, u);
        double b = DoubleArrays.dot(u, v);
        double c = DoubleArrays.dot(v, v);
        double d = DoubleArrays.dot(u, w);
        double e = DoubleArrays.dot(v, w);
        double D = a*c - b*b;
        double sD = D;
        double tD = D;

        double sN, tN;

        if (D < 1e-9) {
            sN = 0.0;
            sD = 1.0;
            tN = e;
            tD = c;
        } else {
            sN = (b*e - c*d);
            tN = (a*e - b*d);
            if (sN < 0) {
                sN = 0;
                tN = e;
                tD = c;
            } else if (sN > sD) {
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }

        if (tN < 0) {
            tN = 0;
            if (-d < 0) {
                sN = 0;
            } else if (-d > a) {
                sN = sD;
            } else {
                sN = -d;
                sD = a;
            }
        } else if (tN > tD) {
            tN = tD;
            if ((-d + b) < 0.0) {
                sN = 0;
            } else if ((-d + b) > a) {
                sN = sD;
            } else {
                sN = (-d + b);
                sD = a;
            }
        }

        double sc = (Math.abs(sN) < 1e-9 ? 0.0 : sN / sD);
        double tc = (Math.abs(tN) < 1e-9 ? 0.0 : tN / tD);

        double sum = 0;
        for (int i = 0; i< dim; ++i) {
            double dPi = w[i] + sc * u[i] - tc * v[i];
            sum += dPi*dPi;
        }

        return Math.sqrt(sum);
    }

    // Prototype method, DO NOT DELETE
    // This is the code before everything is inlined.
    //
//    public static float distSegmentSegment(
//        Vec3f s1p0, Vec3f s1p1,
//        Vec3f s2p0, Vec3f s2p1)
//    {
//        // http://geomalgorithms.com/a07-_distance.html
//        Vec3f u = new Vec3f().sub(s1p1, s1p0);
//        Vec3f v = new Vec3f().sub(s2p1, s2p0);
//        Vec3f w = new Vec3f().sub(s1p0, s2p0);
//
//        float a = u.dot(u);
//        float b = u.dot(v);
//        float c = v.dot(v);
//        float d = u.dot(w);
//        float e = v.dot(w);
//        float D = a*c - b*b;
//        float sD = D;
//        float tD = D;
//
//        float sN, tN;
//
//        if (D < 1e-6f) {
//            sN = 0.0f;
//            sD = 1.0f;
//            tN = e;
//            tD = c;
//        } else {
//            sN = (b*e - c*d);
//            tN = (a*e - b*d);
//            if (sN < 0f) {
//                sN = 0f;
//                tN = e;
//                tD = c;
//            } else if (sN > sD) {
//                sN = sD;
//                tN = e + b;
//                tD = c;
//            }
//        }
//
//        if (tN < 0f) {
//            tN = 0f;
//            if (-d < 0f) {
//                sN = 0f;
//            } else if (-d > a) {
//                sN = sD;
//            } else {
//                sN = -d;
//                sD = a;
//            }
//        } else if (tN > tD) {
//            tN = tD;
//            if ((-d + b) < 0.0f) {
//                sN = 0f;
//            } else if ((-d + b) > a) {
//                sN = sD;
//            } else {
//                sN = (-d + b);
//                sD = a;
//            }
//        }
//
//        float sc = (Math.abs(sN) < 1e-6f ? 0.0f : sN / sD);
//        float tc = (Math.abs(tN) < 1e-6f ? 0.0f : tN / tD);
//
//        // dP = w + (sc * u) - (tc * v);
//        float dPx = w.x + sc * u.x - tc * v.x;
//        float dPy = w.y + sc * u.y - tc * v.y;
//        float dPz = w.z + sc * u.z - tc * v.z;
//
//        return (float)Math.sqrt(dPx*dPx + dPy*dPy + dPz*dPz);
//    }

}
