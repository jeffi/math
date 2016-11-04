// Automatically generated from Vec3d.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

/**
 * Created by jeffi on 11/2/16.
 */
public final class Vec3f implements FloatVector, Cloneable {
    private static final int HASH_BASE = Vec3f.class.getName().hashCode();
    private static final long serialVersionUID = -931885299925000419L;

    public float x;
    public float y;
    public float z;

    public Vec3f() {
    }

    public Vec3f(float x, float y, float z) {
        this.x = x;
        this.y = y;
        this.z = z;
    }

    public Vec3f(float v) {
        this.x = v;
        this.y = v;
        this.z = v;
    }

    public Vec3f(Vec3f v) {
        this.x = v.x;
        this.y = v.y;
        this.z = v.z;
    }

    @Override
    public int columns() {
        return 1;
    }

    @Override
    public int rows() {
        return 3;
    }

    @Override
    public int size() {
        return 0;
    }

    public Vec3f copy(Vec3f v) {
        this.x = v.x;
        this.y = v.y;
        this.z = v.z;
        return this;
    }

    public Vec3f setCoeffs(float x, float y, float z) {
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }

    public Vec3f setCoeffs(float v) {
        this.x = v;
        this.y = v;
        this.z = v;
        return this;
    }

    @Override
    public void setCoeff(int i, float value) {
        switch (i) {
        case 0: x = value; return;
        case 1: y = value; return;
        case 2: z = value; return;
        default:
            throw new MatrixIndexOutOfBoundsException(this, i, 0);
        }
    }

    @Override
    public float getCoeff(int i) {
        switch (i) {
        case 0: return x;
        case 1: return y;
        case 2: return z;
        default:
            throw new MatrixIndexOutOfBoundsException(this, i, 0);
        }
    }

    public Vec3f add(Vec3f a, Vec3f b) {
        this.x = a.x + b.x;
        this.y = a.y + b.y;
        this.z = a.z + b.z;
        return this;
    }

    public Vec3f add(Vec3f v) {
        this.x += v.x;
        this.y += v.y;
        this.z += v.z;
        return this;
    }

    public Vec3f sub(Vec3f a, Vec3f b) {
        this.x = a.x - b.x;
        this.y = a.y - b.y;
        this.z = a.z - b.z;
        return this;
    }

    public Vec3f sub(Vec3f a, float x, float y, float z) {
        this.x = a.x - x;
        this.y = a.y - y;
        this.z = a.z - z;
        return this;
    }

    public Vec3f sub(float x, float y, float z) {
        return sub(this, x, y, z);
    }

    public Vec3f scale(float s) {
        this.x *= s;
        this.y *= s;
        this.z *= s;
        return this;
    }

    public Vec3f scale(Vec3f v, float s) {
        this.x = v.x * s;
        this.y = v.y * s;
        this.z = v.z * s;
        return this;
    }

    public Vec3f transform(AffineTransform3f t, Vec3f v) {
        float vx = v.x;
        float vy = v.y;
        float vz = v.z;
        this.x = t.m00*vx + t.m01*vy + t.m02*vz + t.m03;
        this.y = t.m10*vx + t.m11*vy + t.m12*vz + t.m13;
        this.z = t.m20*vx + t.m21*vy + t.m22*vz + t.m23;
        return this;
    }

    public Vec3f transformUsingRotationOnly(AffineTransform3f t, Vec3f v) {
        float vx = v.x;
        float vy = v.y;
        float vz = v.z;
        this.x = t.m00*vx + t.m01*vy + t.m02*vz;
        this.y = t.m10*vx + t.m11*vy + t.m12*vz;
        this.z = t.m20*vx + t.m21*vy + t.m22*vz;
        return this;
    }

    public Vec3f inverseTransformUsingRotationOnly(AffineTransform3f t, Vec3f v) {
        float vx = v.x;
        float vy = v.y;
        float vz = v.z;
        this.x = t.m00*vx + t.m10*vy + t.m20*vz;
        this.y = t.m01*vx + t.m11*vy + t.m21*vz;
        this.z = t.m02*vx + t.m12*vy + t.m22*vz;
        return this;
    }

    public Vec3f transform(AffineTransform3f t) {
        return this.transform(t, this);
    }

    public Vec3f transform(Quaternion4f q, Vec3f v) {
        // This method uses the formula:
        //    v + 2r x (r x v + wv)
        // Which requires 18 muls + 12 adds = 30 ops
        // It is the fastest transform for quaternions.
        //
        // Matrix transforms however are faster requiring:
        // 9 muls + 6 adds = 15 ops
        //
        // If performing many transforms with the same rotation
        // it is thus better to convert to matrix first, then
        // use the matrix to do the transforms

        final float qw = q.w;
        final float qx = q.x;
        final float qy = q.y;
        final float qz = q.z;
        final float vx = v.x;
        final float vy = v.y;
        final float vz = v.z;

        final float cx = qy*vz - qz*vy + qw * vx;
        final float cy = qz*vx - qx*vz + qw * vy;
        final float cz = qx*vy - qy*vx + qw * vz;

        this.x = 2*(qy*cz - qz*cy) + vx;
        this.y = 2*(qz*cx - qx*cz) + vy;
        this.z = 2*(qx*cy - qy*cx) + vz;
        return this;
    }

    /**
     * Computes this = m * v
     *
     * @param m
     * @param v
     */
    public Vec3f transform(Matrix3f m, Vec3f v) {
        float vx = v.x;
        float vy = v.y;
        float vz = v.z;
        this.x = m.m00*vx + m.m01*vy + m.m02*vz;
        this.y = m.m10*vx + m.m11*vy + m.m12*vz;
        this.z = m.m20*vx + m.m21*vy + m.m22*vz;
        return this;
    }

    /**
     * Transforms a vector by the transpose of a matrix
     * If the matrix is a rotation matrix, the result is
     * the inverse of rotation.  The result is stored in
     * {@code this}.
     *
     * @param m the matrix
     * @param v the vector to transform (may be {@code this})
     * @return {@code this}
     */
    @Deprecated
    public Vec3f inverseTransform(Matrix3f m, Vec3f v) {
        float vx = v.x;
        float vy = v.y;
        float vz = v.z;
        this.x = m.m00*vx + m.m10*vy + m.m20*vz;
        this.y = m.m01*vx + m.m11*vy + m.m21*vz;
        this.z = m.m02*vx + m.m12*vy + m.m22*vz;
        return this;
    }


    /**
     * Computes the transformation of the argument by the inverse of the
     * matrix.  This method should only be used when (1) performing a 1-time
     * inverse (otherwise it would be better to invert the transform and
     * reuse the inverse), and (2) when the transform is not a rigid (translation
     * and rotation only), in which case it would be better to use
     * {@link #fastRigidInverseTransform(AffineTransform3f, Vec3f)}.
     *
     * @param t the transform to use
     * @param v the vector to transform
     * @return this
     */
    public Vec3f inverseTransform(AffineTransform3f t, Vec3f v) {
        // AffineTransform3f it = new AffineTransform3f();
        // it.inverse(t);
        // return this.transform(it, v);

        float a11 = t.m11;
        float a21 = t.m21;
        float a22 = t.m22;
        float a12 = t.m12;
        float inv00 = a11*a22 - a12*a21;
        float a01 = t.m01;
        float a02 = t.m02;
        float inv01 = a02*a21 - a01*a22;
        float inv02 = a01*a12 - a02*a11;
        float a00 = t.m00;
        float a10 = t.m10;
        float a20 = t.m20;
        float det = a00*inv00 + a10*inv01 + a20*inv02;

        if (det == 0)
            throw new ArithmeticException("transform not invertible");

        det = 1.0f/det;

        float dx = v.x - t.m03;
        float dy = v.y - t.m13;
        float dz = v.z - t.m23;
        this.x = (inv00*dx + inv01*dy + inv02*dz)*det;
        float inv10 = a12*a20 - a10*a22;
        float inv11 = a00*a22 - a02*a20;
        float inv12 = a02*a10 - a00*a12;
        this.y = (inv10*dx + inv11*dy + inv12*dz)*det;
        float inv20 = a10*a21 - a11*a20;
        float inv21 = a01*a20 - a00*a21;
        float inv22 = a00*a11 - a01*a10;
        this.z = (inv20*dx + inv21*dy + inv22*dz)*det;
        return this;
    }

    /**
     * Computes the inverse transform of the argument and stores it in
     * {@code this}.  This is equivalent to first inverting the transform
     * then transforming the argument.  To use this method, the transform
     * must be a rigid transform (rotation and translation only).  If there
     * are other transforms (e.g., sheer, scale) in the matrix, then use
     * {@link #inverseTransform(AffineTransform3f, Vec3f)} instead, or
     * first invert the matrix, and transform using the inverted matrix.
     *
     * <p>This method should be preferred in most cases to the slower
     * alternatives.</p>
     *
     * @param t the transform
     * @param v the vector to transform
     * @return this
     */
    public Vec3f fastRigidInverseTransform(AffineTransform3f t, Vec3f v) {
        // AffineTransform3f it = new AffineTransform3f();
        // it.fastRigidInverse(t);
        // return this.transform(it, v);

        final float dx = v.x - t.m03;
        final float dy = v.y - t.m13;
        final float dz = v.z - t.m23;
        this.x = t.m00*dx + t.m10*dy + t.m20*dz;
        this.y = t.m01*dx + t.m11*dy + t.m21*dz;
        this.z = t.m02*dx + t.m12*dy + t.m22*dz;
        return this;
    }

    /**
     * Computes {@code fastRigidInverseTransform(t, this)}.
     *
     * @see #fastRigidInverseTransform(AffineTransform3f, Vec3f)
     * @param t the transform
     * @return {@code this}
     */
    public Vec3f fastRigidInverseTransform(AffineTransform3f t) {
        return fastRigidInverseTransform(t, this);
    }

    /**
     * Computes the dot product of this and another vector.
     *
     * @param other the other vector
     * @return the dot product
     */
    public float dot(Vec3f other) {
        return this.x * other.x + this.y * other.y + this.z * other.z;
    }

    public static float dist(Vec3f a, Vec3f b) {
        float dx = a.x - b.x;
        float dy = a.y - b.y;
        float dz = a.z - b.z;
        return (float)Math.sqrt(dx*dx + dy*dy + dz*dz);
    }

    public static float distSquared(Vec3f a, Vec3f b) {
        float dx = a.x - b.x;
        float dy = a.y - b.y;
        float dz = a.z - b.z;
        return dx*dx + dy*dy + dz*dz;
    }

    public Vec3f cross(Vec3f a, Vec3f b) {
        // compute to temps first in case a == this or b == this
        float x = a.y * b.z - a.z * b.y;
        float y = a.z * b.x - a.x * b.z;
        float z = a.x * b.y - a.y * b.x;
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }

    public float squaredNorm() {
        return x*x + y*y + z*z;
    }

    public float norm() {
        return (float)Math.sqrt(squaredNorm());
    }

    /**
     * Returns the L1 norm of the vector (i.e., the sum of absolute values)
     *
     * @return the L1 norm
     */
    public float L1norm() {
        return Math.abs(x) + Math.abs(y) + Math.abs(z);
    }

    /**
     * Returns the L_infinity norm of the vector (i.e., the max of the
     * absolute values)
     *
     * @return the L_infinity norm
     */
    public float LinfNorm() {
        return Math.max(Math.max(Math.abs(x), Math.abs(y)), Math.abs(z));
    }

//    public boolean normalize() {
//        float s = LinfNorm();
//        float sx = this.x / s;
//        float sy = this.y / s;
//        float sz = this.z / s;
//        float d = (float)(float)Math.sqrt(sx*sx + sy*sy + sz*sz);
//        if (d > 0) {
//            this.x = sx / d;
//            this.y = sy / d;
//            this.z = sz / d;
//            return true;
//        }
//        return false;
//    }

    public boolean normalize(Vec3f v) {
        float d = v.squaredNorm();
        if (d > 0 && d != Float.POSITIVE_INFINITY) {
            d = (float)Math.sqrt(d);
            this.x = v.x / d;
            this.y = v.y / d;
            this.z = v.z / d;
            return true;
        }

        float s = v.LinfNorm();
        if (s > 0) {
            // really small value results in 0 or
            // really large value results in Inf.
            // Pre-scale to get values in a range that we
            // can compute an L2 norm.
            float sx = v.x / s;
            float sy = v.y / s;
            float sz = v.z / s;
            d = (float)Math.sqrt(sx*sx + sy*sy + sz*sz);
            this.x = sx / d;
            this.y = sy / d;
            this.z = sz / d;
            return true;
        }

        return false;
    }

    public boolean normalize() {
        return normalize(this);
    }

    public Vec3f negate(Vec3f v) {
        this.x = -v.x;
        this.y = -v.y;
        this.z = -v.z;
        return this;
    }

    public Vec3f negate() {
        x = -x;
        y = -y;
        z = -z;
        return this;
    }

    public Vec3f min(Vec3f a, Vec3f b) {
        this.x = Math.min(a.x, b.x);
        this.y = Math.min(a.y, b.y);
        this.z = Math.min(a.z, b.z);
        return this;
    }

    public Vec3f min(Vec3f b) {
        return this.min(this, b);
    }

    public Vec3f max(Vec3f a, Vec3f b) {
        this.x = Math.max(a.x, b.x);
        this.y = Math.max(a.y, b.y);
        this.z = Math.max(a.z, b.z);
        return this;
    }

    public Vec3f max(Vec3f b) {
        return this.max(this, b);
    }

    /**
     * Tests if every component of {@code this} is greater
     * than the corresponding component of {@code that}.
     *
     * @param that vector for comparison
     * @return true if every component of {@code this} is greater
     * than the corresponding component of {@code that}
     */
    public boolean allGreaterThan(Vec3f that) {
        return this.x > that.x
               && this.y > that.y
               && this.z > that.z;
    }

    public boolean allGreaterEqual(Vec3f that) {
        return this.x >= that.x
               && this.y >= that.y
               && this.z >= that.z;
    }

    public boolean allLessThan(Vec3f that) {
        return this.x < that.x
               && this.y < that.y
               && this.z < that.z;
    }

    public boolean allLessEqual(Vec3f that) {
        return this.x <= that.x
               && this.y <= that.y
               && this.z <= that.z;
    }


    public float maxCoeff() {
        return Math.max(Math.max(x, y), z);
    }

    public float maxAbsCoeff() {
        return Math.max(Math.max(Math.abs(x), Math.abs(y)), Math.abs(z));
    }

    public Vec3f subScale(Vec3f a, Vec3f b, float s) {
        this.x = a.x - b.x * s;
        this.y = a.y - b.y * s;
        this.z = a.z - b.z * s;
        return this;
    }

    public Vec3f column(Matrix3f m, int col) {
        switch (col) {
        case 0:
            this.x = m.m00;
            this.y = m.m10;
            this.z = m.m20;
            break;
        case 1:
            this.x = m.m01;
            this.y = m.m11;
            this.z = m.m21;
            break;
        case 2:
            this.x = m.m02;
            this.y = m.m12;
            this.z = m.m22;
            break;
        default:
            throw new IndexOutOfBoundsException();
        }
        return this;
    }

    public Vec3f column0(AffineTransform3f m) {
        this.x = m.m00;
        this.y = m.m10;
        this.z = m.m20;
        return this;
    }

    public Vec3f column1(AffineTransform3f m) {
        this.x = m.m01;
        this.y = m.m11;
        this.z = m.m21;
        return this;
    }

    public Vec3f column2(AffineTransform3f m) {
        this.x = m.m02;
        this.y = m.m12;
        this.z = m.m22;
        return this;
    }

    public Vec3f translation(AffineTransform3f t) {
        x = t.m03;
        y = t.m13;
        z = t.m23;
        return this;
    }

    public Vec3f column(AffineTransform3f m, int col) {
        switch (col) {
        case 0:
            return column0(m);
        case 1:
            return column1(m);
        case 2:
            return column2(m);
        case 3:
            return translation(m);
        default:
            throw new IndexOutOfBoundsException();
        }
    }

    public Vec3f addScale(Vec3f a, Vec3f b, float s) {
        this.x = a.x + b.x * s;
        this.y = a.y + b.y * s;
        this.z = a.z + b.z * s;
        return this;
    }

    public Vec3f abs(Vec3f v) {
        this.x = Math.abs(v.x);
        this.y = Math.abs(v.y);
        this.z = Math.abs(v.z);
        return this;
    }

    @Override
    public Vec3f clone() {
        try {
            return (Vec3f)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    @Override
    public int hashCode() {
        // note: 0.0f and -0.0f have different
        // hash codes.  This means that we
        // cannot use Float.compare since it
        // treats 0.0f == -0.0f
        int h = Float.floatToIntBits(x);
        h = h*31 + Float.floatToIntBits(y);
        h = h*31 + Float.floatToIntBits(z);
        return h ^ HASH_BASE;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;
        Vec3f that = (Vec3f)o;
        return MathEx.bitEquals(that.x, x)
               && MathEx.bitEquals(that.y, y)
               && MathEx.bitEquals(that.z, z);
    }

    @Override
    public String toString() {
//        return String.format("(%f, %f, %f)", x, y, z);
        return "(" + x + ", " + y + ", " + z + ")";
    }
}
