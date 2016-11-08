package edu.unc.cs.robotics.math;

/**
 * 2D Geometry methods.
 */
public final class Geom2d {
    private Geom2d() {
        throw new AssertionError("no instances");
    }

    public static double distPointSegmentSquared(
        double px, double py,
        double s0x, double s0y,
        double s1x, double s1y)
    {
        // Vector v = s1 - s0;
        double vx = s1x - s0x;
        double vy = s1y - s0y;
        // Vector w = p - s0;
        double wx = px - s0x;
        double wy = py - s0y;
        // c1 = dot(w, v)
        double c1 = wx * vx + wy * vy;
        if (c1 <= 0) {
            // return dist(p, s0)
            return wx*wx + wy*wy;
        }
        // c2 = dot(v,v)
        double c2 = vx*vx + vy*vy;
        if (c2 <= c1) {
            // return dist(p, s1)
            double ux = px - s1x;
            double uy = py - s1y;
            return ux*ux + uy*uy;
        }
        // b = c1 / c2
        // Pb = s0 + v * b
        // return dist(P, Pb)
        double s = c1 / c2;
        double qx = s0x - px + vx * s;
        double qy = s0y - py + vy * s;
        return qx*qx + qy*qy;
    }

    public static double distPointSegment(
        double px, double py,
        double s0x, double s0y,
        double s1x, double s1y)
    {
        return Math.sqrt(distPointSegmentSquared(px, py, s0x, s0y, s1x, s1y));
    }

    public static double distPointSegmentSquared(Vec2d pt, Vec2d s0, Vec2d s1) {
        return distPointSegmentSquared(pt.x, pt.y, s0.x, s0.y, s1.x, s1.y);
    }


    public static double distPointSegment(Vec2d pt, Vec2d s0, Vec2d s1) {
        return distPointSegment(pt.x, pt.y, s0.x, s0.y, s1.x, s1.y);
    }


    public static double distPointSegmentSquared(
        double px, double py,
        double s0x, double s0y,
        double s1x, double s1y,
        Vec2d nearest)
    {
        // Vector v = s1 - s0;
        double vx = s1x - s0x;
        double vy = s1y - s0y;
        // Vector w = p - s0;
        double wx = px - s0x;
        double wy = py - s0y;
        // c1 = dot(w, v)
        double c1 = wx * vx + wy * vy;
        if (c1 <= 0) {
            // return dist(p, s0)
            if (nearest != null) {
                nearest.x = s0x;
                nearest.y = s0y;
            }
            return wx*wx + wy*wy;
        }
        // c2 = dot(v,v)
        double c2 = vx*vx + vy*vy;
        if (c2 <= c1) {
            // return dist(p, s1)
            double ux = px - s1x;
            double uy = py - s1y;
            if (nearest != null) {
                nearest.x = s1x;
                nearest.y = s1y;
            }
            return ux*ux + uy*uy;
        }
        // b = c1 / c2
        // Pb = s0 + v * b
        // return dist(P, Pb)
        double b = c1 / c2;
        final double nearestX = s0x + vx*b;
        final double nearestY = s0y + vy*b;
        if (nearest != null) {
            nearest.x = nearestX;
            nearest.y = nearestY;
        }
        double qx = nearestX - px;
        double qy = nearestY - py;
        return qx*qx + qy*qy;
    }

    public static double distPointSegment(
        double px, double py,
        double s0x, double s0y,
        double s1x, double s1y,
        Vec2d nearest)
    {
        return Math.sqrt(distPointSegmentSquared(px, py, s0x, s0y, s1x, s1y, nearest));
    }

    public static double distPointSegmentSquared(Vec2d pt, Vec2d s0, Vec2d s1, Vec2d nearest) {
        return distPointSegmentSquared(pt.x, pt.y, s0.x, s0.y, s1.x, s1.y, nearest);
    }

    public static double distPointSegment(Vec2d pt, Vec2d s0, Vec2d s1, Vec2d nearest) {
        return distPointSegment(pt.x, pt.y, s0.x, s0.y, s1.x, s1.y, nearest);
    }

}
