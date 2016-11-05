// Automatically generated from Geom3d.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

/**
 * A collection of geometric computational primitives.
 */
public final class Geom3f {
    private Geom3f() {
        throw new AssertionError("no instances");
    }

    public static float distPointSegmentSquared(Vec3f pt, Vec3f s0, Vec3f s1) {
        return distPointSegmentSquared(pt, s0, s1, null);
    }

    public static float distPointSegmentSquared(Vec3f pt, Vec3f s0, Vec3f s1, Vec3f nearest) {
        // implementation is inlined/optimized from the following:
        //
//        Vec3f v = new Vec3f().sub(s1, s0);
//        Vec3f w = new Vec3f().sub(pt, s0);
//        float c1 = w.dot(v);
//        if (c1 <= 0)
//            return Vec3f.distSquared(pt, s0);
//        float c2 = v.dot(v);
//        if (c2 <= c1)
//            return Vec3f.distSquared(pt, s1);
//        float b = c1 / c2;
//        if (nearest == null)
//            nearest = new Vec3f();
//
//        nearest.mul(v, b).add(s0);
//        return Vec3f.distSquared(pt, nearest);

        float vx = s1.x - s0.x;
        float vy = s1.y - s0.y;
        float vz = s1.z - s0.z;
        float c1 = (pt.x - s0.x)*vx + (pt.y - s0.y)*vy + (pt.z - s0.z)*vz;
        float c2;
        float nx, ny, nz;

        if (c1 <= 0) {
            nx = s0.x;
            ny = s0.y;
            nz = s0.z;
        } else if ((c2 = vx*vx + vy*vy + vz*vz) <= c1) {
            nx = s1.x;
            ny = s1.y;
            nz = s1.z;
        } else {
            float b = c1/c2;
            nx = vx*b + s0.x;
            ny = vy*b + s0.y;
            nz = vz*b + s0.z;
        }

        float dx = pt.x - nx;
        float dy = pt.y - ny;
        float dz = pt.z - nz;
        if (nearest != null) {
            nearest.x = nx;
            nearest.y = ny;
            nearest.z = nz;
        }
        return dx*dx + dy*dy + dz*dz;

    }

    public static float distSquaredSegmentSegment(
        Vec3f s1p0, Vec3f s1p1,
        Vec3f s2p0, Vec3f s2p1)
    {
        return distSquaredSegmentSegment(
            s1p0, s1p1,
            s2p0, s2p1,
            null, null);
    }

    public static float distSquaredSegmentSegment(
        Vec3f s1p0, Vec3f s1p1,
        Vec3f s2p0, Vec3f s2p1,
        Vec3f c1, Vec3f c2)
    {
        final float ux = s1p1.x - s1p0.x;
        final float uy = s1p1.y - s1p0.y;
        final float uz = s1p1.z - s1p0.z;
        final float vx = s2p1.x - s2p0.x;
        final float vy = s2p1.y - s2p0.y;
        final float vz = s2p1.z - s2p0.z;
        final float wx = s1p0.x - s2p0.x;
        final float wy = s1p0.y - s2p0.y;
        final float wz = s1p0.z - s2p0.z;
        final float a = ux*ux + uy*uy + uz*uz;
        final float b = ux*vx + uy*vy + uz*vz;
        final float c = vx*vx + vy*vy + vz*vz;
        final float d = -(ux*wx + uy*wy + uz*wz);
        final float e = vx*wx + vy*wy + vz*wz;
        final float D = a*c - b*b;
        float sD = D;
        float tD = D;

        float sN, tN;

        if (D < 1e-9f) {
            sN = 0.0f;
            sD = 1.0f;
            tN = e;
            tD = c;
        } else {
            sN = (b*e + c*d);
            tN = (a*e + b*d);
            if (sN < 0.0f) {
                sN = 0.0f;
                tN = e;
                tD = c;
            } else if (sN > sD) {
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }

        if (tN < 0.0f) {
            tN = 0.0f;
            if (d < 0.0f) {
                sN = 0.0f;
            } else if (d > a) {
                sN = sD;
            } else {
                sN = d;
                sD = a;
            }
        } else if (tN > tD) {
            tN = tD;
            final float b_d = b + d;
            if (b_d < 0.0f) {
                sN = 0.0f;
            } else if (b_d > a) {
                sN = sD;
            } else {
                sN = b_d;
                sD = a;
            }
        }

        final float sc = (Math.abs(sN) < 1e-9f ? 0.0f : sN/sD);
        final float tc = (Math.abs(tN) < 1e-9f ? 0.0f : tN/tD);

        if (c1 != null && c2 != null) {
            c1.x = s1p1.x + sc * ux;
            c1.y = s1p1.y + sc * uy;
            c1.z = s1p1.z + sc * uz;

            c2.x = s2p1.x + tc * vx;
            c2.y = s2p1.y + tc * vy;
            c2.z = s2p1.z + tc * vz;

            return Vec3f.distSquared(c1, c2);
        } else {
            // if we're not recording the nearest points, we
            // can shortcut w = s1p0 - s2p0 into the distance

            // dP = w + (sc * u) - (tc * v);
            float dPx = wx + sc*ux - tc*vx;
            float dPy = wy + sc*uy - tc*vy;
            float dPz = wz + sc*uz - tc*vz;

            return dPx*dPx + dPy*dPy + dPz*dPz;
        }
    }

    public static float distSegmentSegment(
        Vec3f s1p0, Vec3f s1p1,
        Vec3f s2p0, Vec3f s2p1)
    {
        return (float)Math.sqrt(distSquaredSegmentSegment(s1p0, s1p1, s2p0, s2p1));
    }
}
