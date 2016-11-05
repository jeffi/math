// Automatically generated from Geom3dTest.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/5/16.
 */
public class Geom3fTest extends TestCase {
    static final float EPSILON = 1e-6f;

    private static void checkDistSegmentSegment(
        float expected, float epsilon,
        float ax, float ay, float az,
        float bx, float by, float bz,
        float cx, float cy, float cz,
        float dx, float dy, float dz)
    {
        // checks all 8 variations of the same pair of points.

        final Vec3f a = new Vec3f(ax, ay, az);
        final Vec3f b = new Vec3f(bx, by, bz);
        final Vec3f c = new Vec3f(cx, cy, cz);
        final Vec3f d = new Vec3f(dx, dy, dz);

        float exp = expected;
        float eps = epsilon;

        assertEquals(exp, Geom3f.distSegmentSegment(a, b, c, d), eps);
        assertEquals(exp, Geom3f.distSegmentSegment(b, a, c, d), eps);
        assertEquals(exp, Geom3f.distSegmentSegment(a, b, d, c), eps);
        assertEquals(exp, Geom3f.distSegmentSegment(b, a, d, c), eps);
        assertEquals(exp, Geom3f.distSegmentSegment(c, d, a, b), eps);
        assertEquals(exp, Geom3f.distSegmentSegment(c, d, b, a), eps);
        assertEquals(exp, Geom3f.distSegmentSegment(d, c, a, b), eps);
        assertEquals(exp, Geom3f.distSegmentSegment(d, c, b, a), eps);
    }

    public void testDistSegmentSegment() throws Exception {


        // co-linear segments
        checkDistSegmentSegment(
            0.5f, EPSILON,
            0, 0, 0,
            0, 0, 1,
            0, 0, 1.5f,
            0, 0, 2.5f);

        // overlapping co-linear
        checkDistSegmentSegment(
            0, EPSILON,
            0, 0, 0,
            0, 0, 1,
            0, 0, 0.5f,
            0, 0, 2.5f);

        // perpendicular, point-segment distance
        checkDistSegmentSegment(
            0.7f, EPSILON,
            0, 0, 0,
            0, 0, 4,
            0, 0.7f, 2,
            0, 3, 2);

        // parallel lines
        checkDistSegmentSegment(
            3, EPSILON,
            0, 0, 0,
            0, 0, 4,
            0, 3, 0,
            0, 3, 4);

        // point-point distance in different variations
        checkDistSegmentSegment(
            (float)Math.sqrt(2), EPSILON,
            0, 0, 0,
            0, 0, 4,
            0, 1, -1,
            0, 3, -3);
    }

    private static void checkDistPointSegmentSquared(
        String message,
        float expected, float epsilon,
        float px, float py, float pz,
        float x0, float y0, float z0,
        float x1, float y1, float z1)
    {
        Vec3f pt = new Vec3f(px, py, pz);
        Vec3f s0 = new Vec3f(x0, y0, z0);
        Vec3f s1 = new Vec3f(x1, y1, z1);

        expected *= expected;

        for (int i=0 ; i<3 ; ++i) {
            float actual0 = Geom3f.distPointSegmentSquared(pt, s0, s1);
            assertEquals(message+" shuffle "+i, expected, actual0, epsilon);
            float actual1 = Geom3f.distPointSegmentSquared(pt, s1, s0);
            assertEquals(message+" shuffle "+i, expected, actual1, epsilon);

            shuffle(pt);
            shuffle(s0);
            shuffle(s1);
        }
    }

    private static void shuffle(Vec3f v) {
        float t = v.x;
        v.x = v.y;
        v.y = v.z;
        v.z = t;
    }

    public void testDistPointSegmentSquared() throws Exception {
        checkDistPointSegmentSquared(
            "degenerate segment == point",
            0, 0,
            1, 1, 1,
            1, 1, 1,
            1, 1, 1);

        checkDistPointSegmentSquared(
            "segment endpoint == point",
            0, 0,
            1, 1, 1,
            1, 1, 1,
            2, 2, 2);

        checkDistPointSegmentSquared(
            "point on segment",
            0, 0,
            1, 1, 1,
            -1, -1, -1,
            2, 2, 2);

        checkDistPointSegmentSquared(
            "point nearest segment endpoint",
            (float)Math.sqrt(3), 1e-10f,
            2, 3, 4,
            3, 4, 5,
            4, 8, 12);

        checkDistPointSegmentSquared(
            "point nearest segment center",
            (float)Math.sqrt(5), 1e-10f,
            0, 1, 2,
            -3, 0, 0,
            3, 0, 0);
    }

    public void testDistPointSegmentSquaredPart2() throws Exception {
        // tests from libccd

        Vec3f P = new Vec3f();
        Vec3f a = new Vec3f();
        Vec3f b = new Vec3f();
        Vec3f w = new Vec3f();
        Vec3f ew = new Vec3f();

        float dist;

        a.setCoeffs(0.0f, 0.0f, 0.0f);
        b.setCoeffs(1.0f, 0.0f, 0.0f);

        // extereme w == a
        float x9 = -1.0f;
        P.setCoeffs(x9, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(1.0f, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        float x8 = -0.5f;
        P.setCoeffs(x8, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.5f * 0.5f, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        float x7 = -0.1f;
        P.setCoeffs(x7, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.1f * 0.1f, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        P.setCoeffs(0.0f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.0f, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        float x6 = -1.0f;
        P.setCoeffs(x6, 1.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(2.0f, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        float x5 = -0.5f;
        P.setCoeffs(x5, 0.5f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.5f, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        float x4 = -0.1f;
        float y1 = -1.0f;
        P.setCoeffs(x4, y1, 2.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(5.01f, dist, EPSILON);
        assertEquals(w, a, EPSILON);


        // extereme w == b
        P.setCoeffs(2.0f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(1.0f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(1.5f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.5f * 0.5f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(1.1f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.1f * 0.1f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(1.0f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.0f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(2.0f, 1.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(2.0f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(1.5f, 0.5f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.5f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        float y = -1.0f;
        P.setCoeffs(1.1f, y, 2.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(5.01f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        // inside segment
        P.setCoeffs(0.5f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.0f, dist, EPSILON);
        assertEquals(w, P, EPSILON);

        P.setCoeffs(0.9f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.0f, dist, EPSILON);
        assertEquals(w, P, EPSILON);

        P.setCoeffs(0.5f, 1.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(1.0f, dist, EPSILON);
        ew.setCoeffs(0.5f, 0.0f, 0.0f);
        assertEquals(w, ew, EPSILON);

        P.setCoeffs(0.5f, 1.0f, 1.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(2.0f, dist, EPSILON);
        ew.setCoeffs(0.5f, 0.0f, 0.0f);
        assertEquals(w, ew, EPSILON);


        float x3 = -0.5f;
        a.setCoeffs(x3, 2.0f, 1.0f);
        b.setCoeffs(1.0f, 1.5f, 0.5f);

        // extreme w == a
        float x2 = -10.0f;
        P.setCoeffs(x2, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(9.5f * 9.5f + 2.0f * 2.0f + 1.0f, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        float x1 = -10.0f;
        P.setCoeffs(x1, 9.2f, 3.4f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(9.5f * 9.5f + 7.2f * 7.2f + 2.4f * 2.4f, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        // extereme w == b
        P.setCoeffs(10.0f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(9.0f * 9.0f + 1.5f * 1.5f + 0.5f * 0.5f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(10.0f, 9.2f, 3.4f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(9.0f * 9.0f + 7.7f * 7.7f + 2.9f * 2.9f, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        // inside ab
        float x = -0.1f;
        a.setCoeffs(x, 1.0f, 1.0f);
        b.setCoeffs(1.0f, 1.0f, 1.0f);
        P.setCoeffs(0.0f, 0.0f, 0.0f);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(2.0f, dist, EPSILON);
        ew.setCoeffs(0.0f, 1.0f, 1.0f);
        assertEquals(w, ew, EPSILON);

    }

    public static void assertEquals(Vec3f e, Vec3f a, float delta) {
        assertEquals(".x", e.x, a.x, delta);
        assertEquals(".y", e.y, a.y, delta);
        assertEquals(".z", e.z, a.z, delta);
    }

    private static float ccdVec3PointSegmentDist2(Vec3f p, Vec3f a, Vec3f b, Vec3f w) {
        return Geom3f.distPointSegmentSquared(p, a, b, w);
    }

}