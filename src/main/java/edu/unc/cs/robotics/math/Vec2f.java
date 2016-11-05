// Automatically generated from Vec2d.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

/**
 * Created by jeffi on 11/3/16.
 */
public final class Vec2f implements FloatVector, Cloneable {
    private static final int HASH_BASE = Vec2f.class.getName().hashCode();
    private static final long serialVersionUID = -2162533668147272143L;

    public float x;
    public float y;

    public Vec2f() {}

    public Vec2f(float x, float y) {
        this.x = x;
        this.y = y;
    }

    public Vec2f(float[] v) {
        x = v[0];
        y = v[1];
    }

    public Vec2f set(float x, float y) {
        this.x = x;
        this.y = y;
        return this;
    }

    @Override
    public int rows() {
        return 2;
    }

    @Override
    public int columns() {
        return 1;
    }

    @Override
    public int size() {
        return 2;
    }

    @Override
    public float getCoeff(int i) {
        switch (i) {
        case 0: return x;
        case 1: return y;
        default:
            throw new MatrixIndexOutOfBoundsException(this, i, 0);
        }
    }

    @Override
    public void setCoeff(int i, float v) {
        switch (i) {
        case 0: x = v; return;
        case 1: y = v; return;
        default:
            throw new MatrixIndexOutOfBoundsException(this, i, 0);
        }
    }

    public Vec2f projectX(Vec3f v) {
        return set(v.y, v.z);
    }

    public Vec2f projectY(Vec3f v) {
        return set(v.x, v.z);
    }

    public Vec2f projectZ(Vec3f v) {
        return set(v.x, v.y);
    }

    public Vec2f sub(Vec2f a, Vec2f b) {
        this.x = a.x - b.x;
        this.y = a.y - b.y;
        return this;
    }

    public Vec2f scale(float s) {
        this.x *= s;
        this.y *= s;
        return this;
    }


    public float dot(Vec2f b) {
        return this.x * b.x + this.y * b.y;
    }

    public float squaredNorm() {
        return x*x + y*y;
    }

    public float LinfNorm() {
        return Math.max(Math.abs(x), Math.abs(y));
    }

    public boolean normalize() {
        float d = squaredNorm();
        if (d > 0 && d != Float.POSITIVE_INFINITY) {
            d = (float)Math.sqrt(d);
            this.x /= d;
            this.y /= d;
            return true;
        }

        float s = LinfNorm();
        if (s > 0) {
            // really small value results in 0 or
            // really large value results in Inf.
            // Pre-scale to get values in a range that we
            // can compute an L2 norm.
            float sx = this.x / s;
            float sy = this.y / s;
            d = (float)Math.sqrt(sx*sx + sy*sy);
            this.x = sx / d;
            this.y = sy / d;
            return true;
        }

        return false;
    }

    /**
     * Perpendicular dot product.
     *
     * @param b
     * @return
     */
    public float perpendicularDot(Vec2f b) {
        return this.x * b.y - this.y * b.x;
    }

    @Override
    public Vec2f clone() {
        try {
            return (Vec2f)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    @Override
    public int hashCode() {
        int h = Float.floatToIntBits(x);
        h = h*31 + Float.floatToIntBits(y);
        return h ^ HASH_BASE;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) return true;
        if (obj == null || obj.getClass() != this.getClass()) return false;
        final Vec2f that = (Vec2f)obj;
        return Float.floatToIntBits(this.x) == Float.floatToIntBits(that.x)
            && Float.floatToIntBits(this.y) == Float.floatToIntBits(that.y);
    }

    @Override
    public String toString() {
        return "(" + x + ", " + y + ")";
    }
}
