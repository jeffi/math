package edu.unc.cs.robotics.math;

/**
 * A collection of geometric computational primitives.
 */
public final class Geom3d {
    private Geom3d() {
        throw new AssertionError("no instances");
    }

    public static double distPointSegmentSquared(Vec3d pt, Vec3d s0, Vec3d s1) {
        return distPointSegmentSquared(pt, s0, s1, null);
    }

    public static double distPointSegmentSquared(Vec3d pt, Vec3d s0, Vec3d s1, Vec3d nearest) {
        // implementation is inlined/optimized from the following:
        //
//        Vec3d v = new Vec3d().sub(s1, s0);
//        Vec3d w = new Vec3d().sub(pt, s0);
//        double c1 = w.dot(v);
//        if (c1 <= 0)
//            return Vec3d.distSquared(pt, s0);
//        double c2 = v.dot(v);
//        if (c2 <= c1)
//            return Vec3d.distSquared(pt, s1);
//        double b = c1 / c2;
//        if (nearest == null)
//            nearest = new Vec3d();
//
//        nearest.mul(v, b).add(s0);
//        return Vec3d.distSquared(pt, nearest);

        double vx = s1.x - s0.x;
        double vy = s1.y - s0.y;
        double vz = s1.z - s0.z;
        double c1 = (pt.x - s0.x)*vx + (pt.y - s0.y)*vy + (pt.z - s0.z)*vz;
        double c2;
        double nx, ny, nz;

        if (c1 <= 0) {
            nx = s0.x;
            ny = s0.y;
            nz = s0.z;
        } else if ((c2 = vx*vx + vy*vy + vz*vz) <= c1) {
            nx = s1.x;
            ny = s1.y;
            nz = s1.z;
        } else {
            double b = c1/c2;
            nx = vx*b + s0.x;
            ny = vy*b + s0.y;
            nz = vz*b + s0.z;
        }

        double dx = pt.x - nx;
        double dy = pt.y - ny;
        double dz = pt.z - nz;
        if (nearest != null) {
            nearest.x = nx;
            nearest.y = ny;
            nearest.z = nz;
        }
        return dx*dx + dy*dy + dz*dz;

    }

    public static double distSquaredSegmentSegment(
        Vec3d s1p0, Vec3d s1p1,
        Vec3d s2p0, Vec3d s2p1)
    {
        return distSquaredSegmentSegment(
            s1p0, s1p1,
            s2p0, s2p1,
            null, null);
    }

    public static double distSquaredSegmentSegment(
        Vec3d s1p0, Vec3d s1p1,
        Vec3d s2p0, Vec3d s2p1,
        Vec3d c1, Vec3d c2)
    {
        final double ux = s1p1.x - s1p0.x;
        final double uy = s1p1.y - s1p0.y;
        final double uz = s1p1.z - s1p0.z;
        final double vx = s2p1.x - s2p0.x;
        final double vy = s2p1.y - s2p0.y;
        final double vz = s2p1.z - s2p0.z;
        final double wx = s1p0.x - s2p0.x;
        final double wy = s1p0.y - s2p0.y;
        final double wz = s1p0.z - s2p0.z;
        final double a = ux*ux + uy*uy + uz*uz;
        final double b = ux*vx + uy*vy + uz*vz;
        final double c = vx*vx + vy*vy + vz*vz;
        final double d = -(ux*wx + uy*wy + uz*wz);
        final double e = vx*wx + vy*wy + vz*wz;
        final double D = a*c - b*b;
        double sD = D;
        double tD = D;

        double sN, tN;

        if (D < 1e-9) {
            sN = 0.0;
            sD = 1.0;
            tN = e;
            tD = c;
        } else {
            sN = (b*e + c*d);
            tN = (a*e + b*d);
            if (sN < 0.0) {
                sN = 0.0;
                tN = e;
                tD = c;
            } else if (sN > sD) {
                sN = sD;
                tN = e + b;
                tD = c;
            }
        }

        if (tN < 0.0) {
            tN = 0.0;
            if (d < 0.0) {
                sN = 0.0;
            } else if (d > a) {
                sN = sD;
            } else {
                sN = d;
                sD = a;
            }
        } else if (tN > tD) {
            tN = tD;
            final double b_d = b + d;
            if (b_d < 0.0) {
                sN = 0.0;
            } else if (b_d > a) {
                sN = sD;
            } else {
                sN = b_d;
                sD = a;
            }
        }

        final double sc = (Math.abs(sN) < 1e-9 ? 0.0 : sN/sD);
        final double tc = (Math.abs(tN) < 1e-9 ? 0.0 : tN/tD);

        if (c1 != null && c2 != null) {
            c1.x = s1p1.x + sc * ux;
            c1.y = s1p1.y + sc * uy;
            c1.z = s1p1.z + sc * uz;

            c2.x = s2p1.x + tc * vx;
            c2.y = s2p1.y + tc * vy;
            c2.z = s2p1.z + tc * vz;

            return Vec3d.distSquared(c1, c2);
        } else {
            // if we're not recording the nearest points, we
            // can shortcut w = s1p0 - s2p0 into the distance

            // dP = w + (sc * u) - (tc * v);
            double dPx = wx + sc*ux - tc*vx;
            double dPy = wy + sc*uy - tc*vy;
            double dPz = wz + sc*uz - tc*vz;

            return dPx*dPx + dPy*dPy + dPz*dPz;
        }
    }

    public static double distSegmentSegment(
        Vec3d s1p0, Vec3d s1p1,
        Vec3d s2p0, Vec3d s2p1)
    {
        return Math.sqrt(distSquaredSegmentSegment(s1p0, s1p1, s2p0, s2p1));
    }
}
