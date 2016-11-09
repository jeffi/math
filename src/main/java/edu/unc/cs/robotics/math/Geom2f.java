// Automatically generated from Geom2d.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

/**
 * 2D Geometry methods.
 */
public final class Geom2f {
    private Geom2f() {
        throw new AssertionError("no instances");
    }

    public static float distPointSegmentSquared(
        float px, float py,
        float s0x, float s0y,
        float s1x, float s1y)
    {
        // Vector v = s1 - s0;
        float vx = s1x - s0x;
        float vy = s1y - s0y;
        // Vector w = p - s0;
        float wx = px - s0x;
        float wy = py - s0y;
        // c1 = dot(w, v)
        float c1 = wx * vx + wy * vy;
        if (c1 <= 0) {
            // return dist(p, s0)
            return wx*wx + wy*wy;
        }
        // c2 = dot(v,v)
        float c2 = vx*vx + vy*vy;
        if (c2 <= c1) {
            // return dist(p, s1)
            float ux = px - s1x;
            float uy = py - s1y;
            return ux*ux + uy*uy;
        }
        // b = c1 / c2
        // Pb = s0 + v * b
        // return dist(P, Pb)
        float s = c1 / c2;
        float qx = s0x - px + vx * s;
        float qy = s0y - py + vy * s;
        return qx*qx + qy*qy;
    }

    public static float distPointSegment(
        float px, float py,
        float s0x, float s0y,
        float s1x, float s1y)
    {
        return (float)Math.sqrt(distPointSegmentSquared(px, py, s0x, s0y, s1x, s1y));
    }

    public static float distPointSegmentSquared(Vec2f pt, Vec2f s0, Vec2f s1) {
        return distPointSegmentSquared(pt.x, pt.y, s0.x, s0.y, s1.x, s1.y);
    }


    public static float distPointSegment(Vec2f pt, Vec2f s0, Vec2f s1) {
        return distPointSegment(pt.x, pt.y, s0.x, s0.y, s1.x, s1.y);
    }


    public static float distPointSegmentSquared(
        float px, float py,
        float s0x, float s0y,
        float s1x, float s1y,
        Vec2f nearest)
    {
        // Vector v = s1 - s0;
        float vx = s1x - s0x;
        float vy = s1y - s0y;
        // Vector w = p - s0;
        float wx = px - s0x;
        float wy = py - s0y;
        // c1 = dot(w, v)
        float c1 = wx * vx + wy * vy;
        if (c1 <= 0) {
            // return dist(p, s0)
            if (nearest != null) {
                nearest.x = s0x;
                nearest.y = s0y;
            }
            return wx*wx + wy*wy;
        }
        // c2 = dot(v,v)
        float c2 = vx*vx + vy*vy;
        if (c2 <= c1) {
            // return dist(p, s1)
            float ux = px - s1x;
            float uy = py - s1y;
            if (nearest != null) {
                nearest.x = s1x;
                nearest.y = s1y;
            }
            return ux*ux + uy*uy;
        }
        // b = c1 / c2
        // Pb = s0 + v * b
        // return dist(P, Pb)
        float b = c1 / c2;
        final float nearestX = s0x + vx*b;
        final float nearestY = s0y + vy*b;
        if (nearest != null) {
            nearest.x = nearestX;
            nearest.y = nearestY;
        }
        float qx = nearestX - px;
        float qy = nearestY - py;
        return qx*qx + qy*qy;
    }

    public static float distPointSegment(
        float px, float py,
        float s0x, float s0y,
        float s1x, float s1y,
        Vec2f nearest)
    {
        return (float)Math.sqrt(distPointSegmentSquared(px, py, s0x, s0y, s1x, s1y, nearest));
    }

    public static float distPointSegmentSquared(Vec2f pt, Vec2f s0, Vec2f s1, Vec2f nearest) {
        return distPointSegmentSquared(pt.x, pt.y, s0.x, s0.y, s1.x, s1.y, nearest);
    }

    public static float distPointSegment(Vec2f pt, Vec2f s0, Vec2f s1, Vec2f nearest) {
        return distPointSegment(pt.x, pt.y, s0.x, s0.y, s1.x, s1.y, nearest);
    }

}
