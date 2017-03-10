// Automatically generated from Quaternion4d.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import java.util.Random;

/**
 * A quaternion.  This class focuses on the their use as
 * a representation for rotations in 3D.
 *
 * <p>This class's API generally follows the formula of having {@code this}
 * be the result.  To efficiently use this class, pre-allocate all the
 * result (and operand) quaternions.</p>
 */
public final class Quaternion4f implements FloatVector, Cloneable {
    private static final int HASH_BASE = Quaternion4f.class.getName().hashCode();
    private static final float MAX_NORM_ERROR = 1e-9f;
    private static final long serialVersionUID = -8521535146520428633L;

    public float w;
    public float x;
    public float y;
    public float z;

    public Quaternion4f() {
        w = 1.0f;
    }

    public Quaternion4f(float w, float x, float y, float z) {
        this.w = w;
        this.x = x;
        this.y = y;
        this.z = z;
    }

    /**
     * Creates a quaternion from an array representation in the order
     * {@code {w,x,y,z}}
     *
     * @param q the quaternion in array form.
     */
    public Quaternion4f(float[] q) {
        this.w = q[0];
        this.x = q[1];
        this.y = q[2];
        this.z = q[3];
    }

    public Quaternion4f setIdentity() {
        this.w = 1.0f;
        this.x = 0.0f;
        this.y = 0.0f;
        this.z = 0.0f;
        return this;
    }

    @Override
    public int rows() {
        return 4;
    }

    @Override
    public int columns() {
        return 1;
    }

    @Override
    public int size() {
        return 4;
    }

    public float getCoeff(int i) {
        switch (i) {
        case 0: return w;
        case 1: return x;
        case 2: return y;
        case 3: return z;
        default:
            throw new MatrixIndexOutOfBoundsException(this, 0, i);
        }
    }

    @Override
    public void setCoeff(int i, float v) {
        switch (i) {
        case 0: w = v; break;
        case 1: x = v; break;
        case 2: y = v; break;
        case 3: z = v; break;
        default:
            throw new MatrixIndexOutOfBoundsException(this, 0, i);
        }
    }

    public Quaternion4f fromRPY(float roll, float pitch, float yaw) {
        final float hallRoll = roll / 2.0f;
        float cx = (float)Math.cos(hallRoll);
        float sx = (float)Math.sin(hallRoll);
        final float halfPitch = pitch / 2.0f;
        float cy = (float)Math.cos(halfPitch);
        float sy = (float)Math.sin(halfPitch);
        final float halfYaw = yaw / 2.0f;
        float cz = (float)Math.cos(halfYaw);
        float sz = (float)Math.sin(halfYaw);

        float t0 = cz * cy;
        float t1 = sz * sy;
        float t2 = cz * sy;
        float t3 = sz * cy;

        this.w = t0 * cx + t1 * sx;
        this.x = t0 * sx - t1 * cx;
        this.y = t3 * sx + t2 * cx;
        this.z = t3 * cx - t2 * sx;
        return this;
    }

    public Quaternion4f fromAxisAngle(float x, float y, float z, float angle) {
        float norm = x*x + y*y + z*z;
        if (norm < MAX_NORM_ERROR) {
            setIdentity();
        } else {
            final float halfAngle = angle / 2.0f;
            float s = ((float)Math.sin(halfAngle) / (float)Math.sqrt(norm));
            this.w = (float)Math.cos(halfAngle);
            this.x = x * s;
            this.y = y * s;
            this.z = z * s;
        }
        return this;
    }

    public Quaternion4f fromRotation(float m00, float m01, float m02,
                                     float m10, float m11, float m12,
                                     float m20, float m21, float m22)
    {
        float trace = m00 + m11 + m22;
        if (trace > 0) {
            float s = (float)Math.sqrt(trace + 1.0f);
            this.w = s * 0.5f;
            this.x = (m21 - m12) * (s = 0.5f / s);
            this.y = (m02 - m20) * s;
            this.z = (m10 - m01) * s;
        } else if (m00 > m11 && m00 > m22) {
            float s = (float)Math.sqrt(m00 - m11 - m22 + 1.0f);
            this.x = s * 0.5f;
            this.w = (m21 - m12) * (s = 0.5f / s);
            this.y = (m10 + m01) * s;
            this.z = (m02 + m20) * s;
        } else if (m11 > m22) {
            float s = (float)Math.sqrt(m11 - m00 - m22 + 1.0f);
            this.y = s * 0.5f;
            this.w = (m02 - m20) * (s = 0.5f / s);
            this.x = (m01 + m10) * s;
            this.z = (m12 + m21) * s;
        } else {
            float s = (float)Math.sqrt(m22 - m00 - m11 + 1.0f);
            this.z = s * 0.5f;
            this.w = (m10 - m01) * (s = 0.5f / s);
            this.x = (m02 + m20) * s;
            this.y = (m12 + m21) * s;
        }

        return this;
    }

    public Quaternion4f fromMatrix(Matrix3f m) {
        return fromRotation(
            m.m00, m.m01, m.m02,
            m.m10, m.m11, m.m12,
            m.m20, m.m21, m.m22);
    }

    public Quaternion4f fromRotation(AffineTransform3f t) {
        return fromRotation(
            t.m00, t.m01, t.m02,
            t.m10, t.m11, t.m12,
            t.m20, t.m21, t.m22);
    }

//    public static void fromMatrix(float[] q, float[] m) {
//        fromMatrix(q,
//            m[0], m[3], m[6],
//            m[1], m[4], m[7],
//            m[2], m[5], m[8]);
//    }
//
//    public static void fromMatrix(
//        float[] q,
//        float m00, float m01, float m02,
//        float m10, float m11, float m12,
//        float m20, float m21, float m22)
//    {
//        float trace = m00 + m11 + m22;
//
//        if (trace > 0) {
//            float root = (float)Math.sqrt(trace + 1);
//            q[0] = root * 0.5f;
//            root = 0.5f / root;
//            q[1] = (m21 - m12) * root;
//            q[2] = (m02 - m20) * root;
//            q[3] = (m10 - m01) * root;
//        } else if (m00 > m11 && m00 > m22) {
//            // i = 0, j = 1, k = 2
//            float root = (float)Math.sqrt(m00 - m11 - m22 + 1);
//            q[1] = 0.5f * root;
//            root = 0.5f / root;
//            q[0] = (m21 - m12) * root;
//            q[2] = (m10 + m01) * root;
//            q[3] = (m20 + m02) * root;
//        } else if (m11 > m22) {
//            // i = 1, j = 2, k = 0
//            float root = (float)Math.sqrt(m11 - m22 - m00 + 1);
//            q[2] = 0.5f * root;
//            root = 0.5f / root;
//            q[0] = (m02 - m20) * root;
//            q[3] = (m21 + m12) * root;
//            q[1] = (m01 + m10) * root;
//        } else {
//            // i = 2, j = 0, k = 1
//            float root = (float)Math.sqrt(m22 - m00 - m11 + 1);
//            q[3] = 0.5f * root;
//            root = 0.5f / root;
//            q[0] = (m10 - m01) * root;
//            q[1] = (m02 + m20) * root;
//            q[2] = (m12 + m21) * root;
//        }
//    }
//
//    public static void fromMatrixAlt(
//        float[] q,
//        float m00, float m01, float m02,
//        float m10, float m11, float m12,
//        float m20, float m21, float m22)
//    {
//        // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
//        // "Alternative Method"
//
//        // This method might have problems because of the use of copySign with 0
//        // Specifically, thr following rotation:
//        //   1  0  0
//        //   0  1  0
//        //   0  0 -1
//
//        // The Math.max is a safeguard against rounding error
//        q[0] = (float)Math.sqrt( Math.max(0, 1 + m00 + m11 + m22) ) * 0.5f;
//
//        q[1] = Math.copySign(
//            (float)Math.sqrt( Math.max(0, 1 + m00 - m11 - m22) ) * 0.5f,
//            m21 - m12);
//        q[2] = Math.copySign(
//            (float)Math.sqrt( Math.max(0, 1 - m00 + m11 - m22) ) * 0.5f,
//            m02 - m20);
//        q[3] = Math.copySign(
//            (float)Math.sqrt( Math.max(0, 1 - m00 - m11 + m22) ) * 0.5f,
//            m10 - m01);
//    }
//
//    public static void fromMatrixAlt2(
//        float[] q,
//        float m00, float m01, float m02,
//        float m10, float m11, float m12,
//        float m20, float m21, float m22)
//    {
//        float root = (float)Math.sqrt(1 + m00 + m11 + m22);
//        float w = root * 0.5f;
//        root = 0.5f / root;
//        float x = (m21 - m12) * root;
//        float y = (m02 - m20) * root;
//        float z = (m10 - m01) * root;
//    }

    public float dot(Quaternion4f q) {
        return this.w * q.w
               + this.x * q.x
               + this.y * q.y
               + this.z * q.z;
    }

    /**
     * Computes this = a * b
     * @param a first quaternion operand
     * @param b second quaternion operand
     * @return this
     */
    public Quaternion4f mul(Quaternion4f a, Quaternion4f b) {
        // compute to temporaries first in case r == a || r == b
        float w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
        float x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
        float y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
        float z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;

        this.w = w;
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }

    public Quaternion4f negate(Quaternion4f q) {
        this.w = -q.w;
        this.x = -q.x;
        this.y = -q.y;
        this.z = -q.z;
        return this;
    }

    public Quaternion4f negate() {
        return negate(this);
    }

    public Quaternion4f conj(Quaternion4f q) {
        this.w =  q.w;
        this.x = -q.x;
        this.y = -q.y;
        this.z = -q.z;
        return this;
    }

    public Quaternion4f conj() {
        return conj(this);
    }

    public Quaternion4f negateConj() {
        this.w = -this.w;
        return this;
    }


    public Quaternion4f random(Random random) {
        float x0 = random.nextFloat();
        float r1 = (float)Math.sqrt(1.0f - x0);
        float r2 = (float)Math.sqrt(x0);
        float t1 = random.nextFloat() * ((float)Math.PI * 2);
        float t2 = random.nextFloat() * ((float)Math.PI * 2);
        float c1 = (float)Math.cos(t1);
        float s1 = (float)Math.sin(t1);
        float c2 = (float)Math.cos(t2);
        float s2 = (float)Math.sin(t2);
        this.w = s1 * r1;
        this.x = c1 * r1;
        this.y = s2 * r2;
        this.z = c2 * r2;
        return this;
    }

    public static float arcLength(Quaternion4f q, Quaternion4f p) {
        final float dot = Math.abs(q.w * p.w + q.x * p.x + q.y * p.y + q.z * p.z);
        return dot <= (1.0f - MAX_NORM_ERROR) ? 2.0f * (float)Math.acos(dot) : 0.0f;
    }

    @Override
    public Quaternion4f clone() {
        try {
            return (Quaternion4f)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    @Override
    public int hashCode() {
        int h = Float.floatToIntBits(w);
        h = h*31 + Float.floatToIntBits(x);
        h = h*31 + Float.floatToIntBits(y);
        h = h*31 + Float.floatToIntBits(z);
        return h ^ HASH_BASE;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || obj.getClass() != this.getClass()) return false;
        final Quaternion4f that = (Quaternion4f)obj;
        return MathEx.bitEquals(this.w, that.w)
            && MathEx.bitEquals(this.x, that.x)
            && MathEx.bitEquals(this.y, that.y)
            && MathEx.bitEquals(this.z, that.z);
    }

    @Override
    public String toString() {
        return w + ", " + x + ", " + y + ", " + z;
    }
}
