package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/5/16.
 */
public class Geom3dTest extends TestCase {
    static final double EPSILON = 1e-9;

    private static void checkDistSegmentSegment(
        double expected, double epsilon,
        double ax, double ay, double az,
        double bx, double by, double bz,
        double cx, double cy, double cz,
        double dx, double dy, double dz)
    {
        // checks all 8 variations of the same pair of points.

        final Vec3d a = new Vec3d(ax, ay, az);
        final Vec3d b = new Vec3d(bx, by, bz);
        final Vec3d c = new Vec3d(cx, cy, cz);
        final Vec3d d = new Vec3d(dx, dy, dz);

        double exp = expected;
        double eps = epsilon;

        assertEquals(exp, Geom3d.distSegmentSegment(a, b, c, d), eps);
        assertEquals(exp, Geom3d.distSegmentSegment(b, a, c, d), eps);
        assertEquals(exp, Geom3d.distSegmentSegment(a, b, d, c), eps);
        assertEquals(exp, Geom3d.distSegmentSegment(b, a, d, c), eps);
        assertEquals(exp, Geom3d.distSegmentSegment(c, d, a, b), eps);
        assertEquals(exp, Geom3d.distSegmentSegment(c, d, b, a), eps);
        assertEquals(exp, Geom3d.distSegmentSegment(d, c, a, b), eps);
        assertEquals(exp, Geom3d.distSegmentSegment(d, c, b, a), eps);
    }

    public void testDistSegmentSegment() throws Exception {


        // co-linear segments
        checkDistSegmentSegment(
            0.5, EPSILON,
            0, 0, 0,
            0, 0, 1,
            0, 0, 1.5,
            0, 0, 2.5);

        // overlapping co-linear
        checkDistSegmentSegment(
            0, EPSILON,
            0, 0, 0,
            0, 0, 1,
            0, 0, 0.5,
            0, 0, 2.5);

        // perpendicular, point-segment distance
        checkDistSegmentSegment(
            0.7, EPSILON,
            0, 0, 0,
            0, 0, 4,
            0, 0.7, 2,
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
            Math.sqrt(2), EPSILON,
            0, 0, 0,
            0, 0, 4,
            0, 1, -1,
            0, 3, -3);
    }

    private static void checkDistPointSegmentSquared(
        String message,
        double expected, double epsilon,
        double px, double py, double pz,
        double x0, double y0, double z0,
        double x1, double y1, double z1)
    {
        Vec3d pt = new Vec3d(px, py, pz);
        Vec3d s0 = new Vec3d(x0, y0, z0);
        Vec3d s1 = new Vec3d(x1, y1, z1);

        expected *= expected;

        for (int i=0 ; i<3 ; ++i) {
            double actual0 = Geom3d.distPointSegmentSquared(pt, s0, s1);
            assertEquals(message+" shuffle "+i, expected, actual0, epsilon);
            double actual1 = Geom3d.distPointSegmentSquared(pt, s1, s0);
            assertEquals(message+" shuffle "+i, expected, actual1, epsilon);

            shuffle(pt);
            shuffle(s0);
            shuffle(s1);
        }
    }

    private static void shuffle(Vec3d v) {
        double t = v.x;
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
            Math.sqrt(3), 1e-10,
            2, 3, 4,
            3, 4, 5,
            4, 8, 12);

        checkDistPointSegmentSquared(
            "point nearest segment center",
            Math.sqrt(5), 1e-10,
            0, 1, 2,
            -3, 0, 0,
            3, 0, 0);
    }

    public void testDistPointSegmentSquaredPart2() throws Exception {
        // tests from libccd

        Vec3d P = new Vec3d();
        Vec3d a = new Vec3d();
        Vec3d b = new Vec3d();
        Vec3d w = new Vec3d();
        Vec3d ew = new Vec3d();

        double dist;

        a.setCoeffs(0.0, 0.0, 0.0);
        b.setCoeffs(1.0, 0.0, 0.0);

        // extereme w == a
        double x9 = -1.0;
        P.setCoeffs(x9, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(1.0, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        double x8 = -0.5;
        P.setCoeffs(x8, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.5 * 0.5, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        double x7 = -0.1;
        P.setCoeffs(x7, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.1 * 0.1, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        P.setCoeffs(0.0, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.0, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        double x6 = -1.0;
        P.setCoeffs(x6, 1.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(2.0, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        double x5 = -0.5;
        P.setCoeffs(x5, 0.5, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.5, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        double x4 = -0.1;
        double y1 = -1.0;
        P.setCoeffs(x4, y1, 2.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(5.01, dist, EPSILON);
        assertEquals(w, a, EPSILON);


        // extereme w == b
        P.setCoeffs(2.0, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(1.0, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(1.5, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.5 * 0.5, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(1.1, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.1 * 0.1, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(1.0, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.0, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(2.0, 1.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(2.0, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(1.5, 0.5, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.5, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        double y = -1.0;
        P.setCoeffs(1.1, y, 2.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(5.01, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        // inside segment
        P.setCoeffs(0.5, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.0, dist, EPSILON);
        assertEquals(w, P, EPSILON);

        P.setCoeffs(0.9, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(0.0, dist, EPSILON);
        assertEquals(w, P, EPSILON);

        P.setCoeffs(0.5, 1.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(1.0, dist, EPSILON);
        ew.setCoeffs(0.5, 0.0, 0.0);
        assertEquals(w, ew, EPSILON);

        P.setCoeffs(0.5, 1.0, 1.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(2.0, dist, EPSILON);
        ew.setCoeffs(0.5, 0.0, 0.0);
        assertEquals(w, ew, EPSILON);


        double x3 = -0.5;
        a.setCoeffs(x3, 2.0, 1.0);
        b.setCoeffs(1.0, 1.5, 0.5);

        // extreme w == a
        double x2 = -10.0;
        P.setCoeffs(x2, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(9.5 * 9.5 + 2.0 * 2.0 + 1.0, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        double x1 = -10.0;
        P.setCoeffs(x1, 9.2, 3.4);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(9.5 * 9.5 + 7.2 * 7.2 + 2.4 * 2.4, dist, EPSILON);
        assertEquals(w, a, EPSILON);

        // extereme w == b
        P.setCoeffs(10.0, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(9.0 * 9.0 + 1.5 * 1.5 + 0.5 * 0.5, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        P.setCoeffs(10.0, 9.2, 3.4);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(9.0 * 9.0 + 7.7 * 7.7 + 2.9 * 2.9, dist, EPSILON);
        assertEquals(w, b, EPSILON);

        // inside ab
        double x = -0.1;
        a.setCoeffs(x, 1.0, 1.0);
        b.setCoeffs(1.0, 1.0, 1.0);
        P.setCoeffs(0.0, 0.0, 0.0);
        dist = ccdVec3PointSegmentDist2(P, a, b, w);
        assertEquals(2.0, dist, EPSILON);
        ew.setCoeffs(0.0, 1.0, 1.0);
        assertEquals(w, ew, EPSILON);

    }

    public static void assertEquals(Vec3d e, Vec3d a, double delta) {
        assertEquals(".x", e.x, a.x, delta);
        assertEquals(".y", e.y, a.y, delta);
        assertEquals(".z", e.z, a.z, delta);
    }

    private static double ccdVec3PointSegmentDist2(Vec3d p, Vec3d a, Vec3d b, Vec3d w) {
        return Geom3d.distPointSegmentSquared(p, a, b, w);
    }

}