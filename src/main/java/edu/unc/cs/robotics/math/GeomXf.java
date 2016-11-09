// Automatically generated from GeomXd.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

/**
 * Geometry methods that work for arbitrary dimensions (or at least 2 and 3),
 * and operate on float[].  This should be deprecated.
 */
public final class GeomXf {

    public static final float EPSILON = 1e-6f;

    private GeomXf() {
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
    public static float cpaTime(float[] p1, float[] v1,
                                 float[] p2, float[] v2)
    {
        // t_cpa = -w_0 * (u - v) / |u - v|^2
        // w_0 = P_0 - Q_0

        float dv2 = 0;
        float dot = 0;

        for (int i = p1.length ; --i >= 0 ; ) {
            final float dvi = v1[i] - v2[i];
            dv2 += dvi*dvi;
            final float w0i = p2[i] - p1[i];
            dot += w0i*dvi;
        }

        return dv2 < EPSILON ? 0 : dot/dv2;
    }

    static class DistSegmentSegment {
        final int dim;
        float[] u;
        float[] v;
        float[] w;

        DistSegmentSegment(int dim) {
            this.dim = dim;
            u = new float[dim];
            v = new float[dim];
            w = new float[dim];
        }

        float compute(float[] s1p0, float[] s1p1,
                       float[] s2p0, float[] s2p1)
        {
            FloatArrays.sub(u, s1p1, s1p0);
            FloatArrays.sub(v, s2p1, s2p0);
            FloatArrays.sub(w, s1p0, s2p0);

            float a = FloatArrays.dot(u, u);
            float b = FloatArrays.dot(u, v);
            float c = FloatArrays.dot(v, v);
            float d = FloatArrays.dot(u, w);
            float e = FloatArrays.dot(v, w);
            float D = a*c - b*b;
            float sD = D;
            float tD = D;

            float sN, tN;

            if (D < 1e-9f) {
                sN = 0.0f;
                sD = 1.0f;
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
                if ((-d + b) < 0.0f) {
                    sN = 0;
                } else if ((-d + b) > a) {
                    sN = sD;
                } else {
                    sN = (-d + b);
                    sD = a;
                }
            }

            float sc = (Math.abs(sN) < 1e-9f ? 0.0f : sN / sD);
            float tc = (Math.abs(tN) < 1e-9f ? 0.0f : tN / tD);

            float sum = 0;
            for (int i=0 ; i<dim ; ++i) {
                float dPi = w[i] + sc * u[i] - tc * v[i];
                sum += dPi*dPi;
            }

            return (float)Math.sqrt(sum);
        }
    }

    public static float distPointSegment(float[] pt, float[] s0, float[] s1) {
        final int dim = pt.length;

        float[] v = new float[dim];
        FloatArrays.sub(v, s1, s0);
        float[] w = new float[dim];
        FloatArrays.sub(w, pt, s0);
        float c1 = FloatArrays.dot(w, v);
        float c2;
        if (c1 <= 0) {
            return FloatArrays.length(w);
        } else if ((c2 = FloatArrays.dot(v, v)) <= c1) {
            return FloatArrays.dist(pt, s1);
        } else {
            float b = c1 / c2;
            float sum = 0;
            for (int i = dim; --i >= 0; ) {
                float di = s0[i] + v[i] * b - pt[i];
                sum += di * di;
            }
            return (float)Math.sqrt(sum);
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
    public static float distPointSegmentSquared(float[] pt, float[] s0, float[] s1, float[] nearest) {
        final int dim = pt.length;
        float c1 = 0;
        float c2 = 0;
        // float ww = 0;
        for (int i=dim ; --i >= 0 ; ) {
            float vi = s1[i] - s0[i];
            float wi = pt[i] - s0[i];
            c1 += vi * wi;
            c2 += vi * vi;
            // ww += wi * wi;
            // borrow nearest for temporary storage.
            nearest[i] = vi;
        }
        if (c1 <= 0) {
            System.arraycopy(s0, 0, nearest, 0, dim);
            return FloatArrays.distSquared(pt, s0);
        }
        if (c2 <= c1) {
            System.arraycopy(s1, 0, nearest, 0, dim);
            return FloatArrays.distSquared(pt, s1);
        }
        float b = c1 / c2;
        float sum = 0;
        for (int i=dim ; --i >= 0 ; ) {
            nearest[i] = s0[i] + nearest[i] * b;
            float di = nearest[i] - pt[i];
            sum += di*di;
        }
        return sum;
    }

    public static float distSegmentSegment(
        float[] s1p0, float[] s1p1,
        float[] s2p0, float[] s2p1)
    {
        final int dim = s1p0.length;
        float[] u = new float[dim];
        float[] v = new float[dim];
        float[] w = new float[dim];
        FloatArrays.sub(u, s1p1, s1p0);
        FloatArrays.sub(v, s2p1, s2p0);
        FloatArrays.sub(w, s1p0, s2p0);

        float a = FloatArrays.dot(u, u);
        float b = FloatArrays.dot(u, v);
        float c = FloatArrays.dot(v, v);
        float d = FloatArrays.dot(u, w);
        float e = FloatArrays.dot(v, w);
        float D = a*c - b*b;
        float sD = D;
        float tD = D;

        float sN, tN;

        if (D < 1e-9f) {
            sN = 0.0f;
            sD = 1.0f;
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
            if ((-d + b) < 0.0f) {
                sN = 0;
            } else if ((-d + b) > a) {
                sN = sD;
            } else {
                sN = (-d + b);
                sD = a;
            }
        }

        float sc = (Math.abs(sN) < 1e-9f ? 0.0f : sN / sD);
        float tc = (Math.abs(tN) < 1e-9f ? 0.0f : tN / tD);

        float sum = 0;
        for (int i = 0; i< dim; ++i) {
            float dPi = w[i] + sc * u[i] - tc * v[i];
            sum += dPi*dPi;
        }

        return (float)Math.sqrt(sum);
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
//        if (D < 1e-6ff) {
//            sN = 0.0ff;
//            sD = 1.0ff;
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
//            if ((-d + b) < 0.0ff) {
//                sN = 0f;
//            } else if ((-d + b) > a) {
//                sN = sD;
//            } else {
//                sN = (-d + b);
//                sD = a;
//            }
//        }
//
//        float sc = (Math.abs(sN) < 1e-6ff ? 0.0ff : sN / sD);
//        float tc = (Math.abs(tN) < 1e-6ff ? 0.0ff : tN / tD);
//
//        // dP = w + (sc * u) - (tc * v);
//        float dPx = w.x + sc * u.x - tc * v.x;
//        float dPy = w.y + sc * u.y - tc * v.y;
//        float dPz = w.z + sc * u.z - tc * v.z;
//
//        return (float)(float)Math.sqrt(dPx*dPx + dPy*dPy + dPz*dPz);
//    }

}
