package edu.unc.cs.robotics.math;

/**
 * Created by jeffi on 11/3/16.
 */
public final class Vec2d implements DoubleVector, Cloneable {
    private static final int HASH_BASE = Vec2d.class.getName().hashCode();
    private static final long serialVersionUID = -2162533668147272143L;

    public double x;
    public double y;

    public Vec2d() {}

    public Vec2d(double x, double y) {
        this.x = x;
        this.y = y;
    }

    public Vec2d set(double x, double y) {
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
    public double getCoeff(int i) {
        switch (i) {
        case 0: return x;
        case 1: return y;
        default:
            throw new MatrixIndexOutOfBoundsException(this, i, 0);
        }
    }

    @Override
    public void setCoeff(int i, double v) {
        switch (i) {
        case 0: x = v; return;
        case 1: y = v; return;
        default:
            throw new MatrixIndexOutOfBoundsException(this, i, 0);
        }
    }

    public Vec2d projectX(Vec3d v) {
        return set(v.y, v.z);
    }

    public Vec2d projectY(Vec3d v) {
        return set(v.x, v.z);
    }

    public Vec2d projectZ(Vec3d v) {
        return set(v.x, v.y);
    }

    public Vec2d sub(Vec2d a, Vec2d b) {
        this.x = a.x - b.x;
        this.y = a.y - b.y;
        return this;
    }

    public Vec2d scale(double s) {
        this.x *= s;
        this.y *= s;
        return this;
    }


    public double dot(Vec2d b) {
        return this.x * b.x + this.y * b.y;
    }

    public double squaredNorm() {
        return x*x + y*y;
    }

    public double LinfNorm() {
        return Math.max(Math.abs(x), Math.abs(y));
    }

    public boolean normalize() {
        double d = squaredNorm();
        if (d > 0 && d != Double.POSITIVE_INFINITY) {
            d = Math.sqrt(d);
            this.x /= d;
            this.y /= d;
            return true;
        }

        double s = LinfNorm();
        if (s > 0) {
            // really small value results in 0 or
            // really large value results in Inf.
            // Pre-scale to get values in a range that we
            // can compute an L2 norm.
            double sx = this.x / s;
            double sy = this.y / s;
            d = Math.sqrt(sx*sx + sy*sy);
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
    public double perpendicularDot(Vec2d b) {
        return this.x * b.y - this.y * b.x;
    }

    @Override
    public Vec2d clone() {
        try {
            return (Vec2d)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    @Override
    public int hashCode() {
        long h = Double.doubleToLongBits(x);
        h = h*31 + Double.doubleToLongBits(y);
        return (int)h ^ (int)(h >>> 32) ^ HASH_BASE;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) return true;
        if (obj == null || obj.getClass() != this.getClass()) return false;
        final Vec2d that = (Vec2d)obj;
        return Double.doubleToLongBits(this.x) == Double.doubleToLongBits(that.x)
            && Double.doubleToLongBits(this.y) == Double.doubleToLongBits(that.y);
    }

    @Override
    public String toString() {
        return "(" + x + ", " + y + ")";
    }
}
