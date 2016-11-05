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
public final class Quaternion4d implements DoubleVector, Cloneable {
    private static final int HASH_BASE = Quaternion4d.class.getName().hashCode();
    private static final double MAX_NORM_ERROR = 1e-9;
    private static final long serialVersionUID = -8521535146520428633L;

    public double w;
    public double x;
    public double y;
    public double z;

    public Quaternion4d() {
        w = 1.0;
    }

    public Quaternion4d(double w, double x, double y, double z) {
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
    public Quaternion4d(double[] q) {
        this.w = q[0];
        this.x = q[1];
        this.y = q[2];
        this.z = q[3];
    }

    public Quaternion4d setIdentity() {
        this.w = 1.0;
        this.x = 0.0;
        this.y = 0.0;
        this.z = 0.0;
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

    public double getCoeff(int i) {
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
    public void setCoeff(int i, double v) {
        switch (i) {
        case 0: w = v; break;
        case 1: x = v; break;
        case 2: y = v; break;
        case 3: z = v; break;
        default:
            throw new MatrixIndexOutOfBoundsException(this, 0, i);
        }
    }

    public Quaternion4d fromRPY(double roll, double pitch, double yaw) {
        final double hallRoll = roll / 2.0;
        double cx = Math.cos(hallRoll);
        double sx = Math.sin(hallRoll);
        final double halfPitch = pitch / 2.0;
        double cy = Math.cos(halfPitch);
        double sy = Math.sin(halfPitch);
        final double halfYaw = yaw / 2.0;
        double cz = Math.cos(halfYaw);
        double sz = Math.sin(halfYaw);

        double t0 = cz * cy;
        double t1 = sz * sy;
        double t2 = cz * sy;
        double t3 = sz * cy;

        this.w = t0 * cx + t1 * sx;
        this.x = t0 * sx - t1 * cx;
        this.y = t3 * sx + t2 * cx;
        this.z = t3 * cx - t2 * sx;
        return this;
    }

    public Quaternion4d fromAxisAngle(double x, double y, double z, double angle) {
        double norm = x*x + y*y + z*z;
        if (norm < MAX_NORM_ERROR) {
            setIdentity();
        } else {
            final double halfAngle = angle / 2.0;
            double s = (Math.sin(halfAngle) / Math.sqrt(norm));
            this.w = Math.cos(halfAngle);
            this.x = x * s;
            this.y = y * s;
            this.z = z * s;
        }
        return this;
    }

    public Quaternion4d fromRotation(double m00, double m01, double m02,
                                     double m10, double m11, double m12,
                                     double m20, double m21, double m22)
    {
        double trace = m00 + m11 + m22;
        if (trace > 0) {
            double s = Math.sqrt(trace + 1.0);
            this.w = s * 0.5;
            this.x = (m21 - m12) * (s = 0.5 / s);
            this.y = (m02 - m20) * s;
            this.z = (m10 - m01) * s;
        } else if (m00 > m11 && m00 > m22) {
            double s = Math.sqrt(m00 - m11 - m22 + 1.0);
            this.x = s * 0.5;
            this.w = (m21 - m12) * (s = 0.5 / s);
            this.y = (m10 + m01) * s;
            this.z = (m02 + m20) * s;
        } else if (m11 > m22) {
            double s = Math.sqrt(m11 - m00 - m22 + 1.0);
            this.y = s * 0.5;
            this.w = (m02 - m20) * (s = 0.5 / s);
            this.x = (m01 + m10) * s;
            this.z = (m12 + m21) * s;
        } else {
            double s = Math.sqrt(m22 - m00 - m11 + 1.0);
            this.z = s * 0.5;
            this.w = (m10 - m01) * (s = 0.5 / s);
            this.x = (m02 + m20) * s;
            this.y = (m12 + m21) * s;
        }

        return this;
    }

    public Quaternion4d fromMatrix(Matrix3d m) {
        return fromRotation(
            m.m00, m.m01, m.m02,
            m.m10, m.m11, m.m12,
            m.m20, m.m21, m.m22);
    }

    public Quaternion4d fromRotation(AffineTransform3d t) {
        return fromRotation(
            t.m00, t.m01, t.m02,
            t.m10, t.m11, t.m12,
            t.m20, t.m21, t.m22);
    }

//    public static void fromMatrix(double[] q, double[] m) {
//        fromMatrix(q,
//            m[0], m[3], m[6],
//            m[1], m[4], m[7],
//            m[2], m[5], m[8]);
//    }
//
//    public static void fromMatrix(
//        double[] q,
//        double m00, double m01, double m02,
//        double m10, double m11, double m12,
//        double m20, double m21, double m22)
//    {
//        double trace = m00 + m11 + m22;
//
//        if (trace > 0) {
//            double root = Math.sqrt(trace + 1);
//            q[0] = root * 0.5;
//            root = 0.5 / root;
//            q[1] = (m21 - m12) * root;
//            q[2] = (m02 - m20) * root;
//            q[3] = (m10 - m01) * root;
//        } else if (m00 > m11 && m00 > m22) {
//            // i = 0, j = 1, k = 2
//            double root = Math.sqrt(m00 - m11 - m22 + 1);
//            q[1] = 0.5 * root;
//            root = 0.5 / root;
//            q[0] = (m21 - m12) * root;
//            q[2] = (m10 + m01) * root;
//            q[3] = (m20 + m02) * root;
//        } else if (m11 > m22) {
//            // i = 1, j = 2, k = 0
//            double root = Math.sqrt(m11 - m22 - m00 + 1);
//            q[2] = 0.5 * root;
//            root = 0.5 / root;
//            q[0] = (m02 - m20) * root;
//            q[3] = (m21 + m12) * root;
//            q[1] = (m01 + m10) * root;
//        } else {
//            // i = 2, j = 0, k = 1
//            double root = Math.sqrt(m22 - m00 - m11 + 1);
//            q[3] = 0.5 * root;
//            root = 0.5 / root;
//            q[0] = (m10 - m01) * root;
//            q[1] = (m02 + m20) * root;
//            q[2] = (m12 + m21) * root;
//        }
//    }
//
//    public static void fromMatrixAlt(
//        double[] q,
//        double m00, double m01, double m02,
//        double m10, double m11, double m12,
//        double m20, double m21, double m22)
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
//        q[0] = Math.sqrt( Math.max(0, 1 + m00 + m11 + m22) ) * 0.5;
//
//        q[1] = Math.copySign(
//            Math.sqrt( Math.max(0, 1 + m00 - m11 - m22) ) * 0.5,
//            m21 - m12);
//        q[2] = Math.copySign(
//            Math.sqrt( Math.max(0, 1 - m00 + m11 - m22) ) * 0.5,
//            m02 - m20);
//        q[3] = Math.copySign(
//            Math.sqrt( Math.max(0, 1 - m00 - m11 + m22) ) * 0.5,
//            m10 - m01);
//    }
//
//    public static void fromMatrixAlt2(
//        double[] q,
//        double m00, double m01, double m02,
//        double m10, double m11, double m12,
//        double m20, double m21, double m22)
//    {
//        double root = Math.sqrt(1 + m00 + m11 + m22);
//        double w = root * 0.5;
//        root = 0.5 / root;
//        double x = (m21 - m12) * root;
//        double y = (m02 - m20) * root;
//        double z = (m10 - m01) * root;
//    }

    public double dot(Quaternion4d q) {
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
    public Quaternion4d mul(Quaternion4d a, Quaternion4d b) {
        // compute to temporaries first in case r == a || r == b
        double w = a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z;
        double x = a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y;
        double y = a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x;
        double z = a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w;

        this.w = w;
        this.x = x;
        this.y = y;
        this.z = z;
        return this;
    }

    public Quaternion4d negate(Quaternion4d q) {
        this.w = -q.w;
        this.x = -q.x;
        this.y = -q.y;
        this.z = -q.z;
        return this;
    }

    public Quaternion4d negate() {
        return negate(this);
    }

    public Quaternion4d conj(Quaternion4d q) {
        this.w =  q.w;
        this.x = -q.x;
        this.y = -q.y;
        this.z = -q.z;
        return this;
    }

    public Quaternion4d conj() {
        return conj(this);
    }

    public Quaternion4d negateConj() {
        this.w = -this.w;
        return this;
    }


    public Quaternion4d random(Random random) {
        double x0 = random.nextDouble();
        double r1 = Math.sqrt(1.0 - x0);
        double r2 = Math.sqrt(x0);
        double t1 = random.nextDouble() * (Math.PI * 2);
        double t2 = random.nextDouble() * (Math.PI * 2);
        double c1 = Math.cos(t1);
        double s1 = Math.sin(t1);
        double c2 = Math.cos(t2);
        double s2 = Math.sin(t2);
        this.w = s1 * r1;
        this.x = c1 * r1;
        this.y = s2 * r2;
        this.z = c2 * r2;
        return this;
    }

    public static double arcLength(Quaternion4d q, Quaternion4d p) {
        final double dot = Math.abs(q.w * p.w + q.x * p.x + q.y * p.y + q.z * p.z);
        return dot > (1.0 - MAX_NORM_ERROR) ? 0.0 : Math.acos(dot);
    }

    @Override
    public Quaternion4d clone() {
        try {
            return (Quaternion4d)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    @Override
    public int hashCode() {
        long h = Double.doubleToLongBits(w);
        h = h*31 + Double.doubleToLongBits(x);
        h = h*31 + Double.doubleToLongBits(y);
        h = h*31 + Double.doubleToLongBits(z);
        return (int)h ^ (int)(h >>> 32) ^ HASH_BASE;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (obj == null || obj.getClass() != this.getClass()) return false;
        final Quaternion4d that = (Quaternion4d)obj;
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
