// Automatically generated from Matrix3d.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import java.util.function.DoubleFunction;

/**
 * A 3x3 matrix, specialized for representing rotations in 3D.
 *
 * <p>This class's API generally follows the formula of having {@code this}
 * be the result.  To efficiently use this class, pre-allocate all the
 * result (and operand) matrices.</p>
 */
public final class Matrix3f implements FloatMatrix, Cloneable {

    private static final int HASH_BASE = Matrix3f.class.getName().hashCode();
    private static final float EPSILON = 1e-6f;

    private static final long serialVersionUID = -2720018621259089452L;

    public float m00;
    public float m01;
    public float m02;
    public float m10;
    public float m11;
    public float m12;
    public float m20;
    public float m21;
    public float m22;

    public Matrix3f() {
        m00 = 1.0f;
        m11 = 1.0f;
        m22 = 1.0f;
    }

    public Matrix3f(Matrix3f m) {
        this.m00 = m.m00;
        this.m01 = m.m01;
        this.m02 = m.m02;
        this.m10 = m.m10;
        this.m11 = m.m11;
        this.m12 = m.m12;
        this.m20 = m.m20;
        this.m21 = m.m21;
        this.m22 = m.m22;
    }

    public Matrix3f(
        float m00, float m01, float m02,
        float m10, float m11, float m12,
        float m20, float m21, float m22)
    {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
    }

    /**
     * Creates a matrix from the coefficients in a column-major array.
     *
     * @param m the column-major ordered matrix
     */
    public Matrix3f(float[] m) {
        m00 = m[0];
        m10 = m[1];
        m20 = m[2];
        m01 = m[3];
        m11 = m[4];
        m21 = m[5];
        m02 = m[6];
        m12 = m[7];
        m22 = m[8];
    }

    /**
     * Sets the matrix to the specified values.  Argument order
     * is in row-order to match visually with the matrix.
     *
     * @param m00
     * @param m01
     * @param m02
     * @param m10
     * @param m11
     * @param m12
     * @param m20
     * @param m21
     * @param m22
     */
    public Matrix3f setCoeffs(
        float m00, float m01, float m02,
        float m10, float m11, float m12,
        float m20, float m21, float m22)
    {
        this.m00 = m00;
        this.m01 = m01;
        this.m02 = m02;
        this.m10 = m10;
        this.m11 = m11;
        this.m12 = m12;
        this.m20 = m20;
        this.m21 = m21;
        this.m22 = m22;
        return this;
    }

    @Override
    public int rows() {
        return 3;
    }

    @Override
    public int columns() {
        return 3;
    }

    @Override
    public int size() {
        return 9;
    }

    public Matrix3f setIdentity() {
        this.m00 = 1.0f;
        this.m01 = 0.0f;
        this.m02 = 0.0f;
        this.m10 = 0.0f;
        this.m11 = 1.0f;
        this.m12 = 0.0f;
        this.m20 = 0.0f;
        this.m21 = 0.0f;
        this.m22 = 1.0f;
        return this;
    }

    public Matrix3f fromTransform(AffineTransform3f t) {
        this.m00 = t.m00;
        this.m10 = t.m10;
        this.m20 = t.m20;
        this.m01 = t.m01;
        this.m11 = t.m11;
        this.m21 = t.m21;
        this.m02 = t.m02;
        this.m12 = t.m12;
        this.m22 = t.m22;
        return this;
    }

    @Override
    public float getCoeff(int r, int c) {
        if (c < 0 || c > 2)
            throw new IndexOutOfBoundsException();

        float v;

        switch (r*3 + c) {
        case 0: v = m00; break;
        case 1: v = m01; break;
        case 2: v = m02; break;
        case 3: v = m10; break;
        case 4: v = m11; break;
        case 5: v = m12; break;
        case 6: v = m20; break;
        case 7: v = m21; break;
        case 8: v = m22; break;
        default:
            throw new IndexOutOfBoundsException();
        }

        return v;
    }

    @Override
    public void setCoeff(int r, int c, float v) {
        if (c < 0 || c > 2)
            throw new IndexOutOfBoundsException();

        switch (r*3+c) {
        case 0: m00 = v; break;
        case 1: m01 = v; break;
        case 2: m02 = v; break;
        case 3: m10 = v; break;
        case 4: m11 = v; break;
        case 5: m12 = v; break;
        case 6: m20 = v; break;
        case 7: m21 = v; break;
        case 8: m22 = v; break;
        default:
            throw new IndexOutOfBoundsException();
        }
    }

    public Matrix3f transpose(Matrix3f a) {
        // temporary for the case r == a
        float tmp;
        this.m00 = a.m00;
        this.m11 = a.m11;
        this.m22 = a.m22;
        tmp = a.m01;
        this.m01 = a.m10;
        this.m10 = tmp;
        tmp = a.m02;
        this.m02 = a.m20;
        this.m20 = tmp;
        tmp = a.m12;
        this.m12 = a.m21;
        this.m21 = tmp;
        return this;
    }

    public Matrix3f transpose() {
        return transpose(this);
    }

    public Matrix3f setEulerRPY(float x, float y, float z) {
        float cx = (float)Math.cos(x),
            sx = (float)Math.sin(x),
            cy = (float)Math.cos(y),
            sy = (float)Math.sin(y),
            cz = (float)Math.cos(z),
            sz = (float)Math.sin(z),
            cc = cx * cz,
            cs = cx * sz,
            sc = sx * cz,
            ss = sx * sz;
        this.m00 = cy * cz;
        this.m10 = cy * sz;
        this.m20 = -sy;
        this.m01 = sy * sc - cs;
        this.m11 = sy * ss + cc;
        this.m21 = cy * sx;
        this.m02 = sy * cc + ss;
        this.m12 = sy * cs - sc;
        this.m22 = cy * cx;
        return this;
    }

    public Matrix3f fromQuaternion(Quaternion4f q) {
        float w2 = q.w*2;
        float x = q.x;
        float y = q.y;
        float z = q.z;
        float wx = w2*x;
        float wy = w2*y;
        float wz = w2*z;
        float x2 = x*2;
        float xx = x2*x;
        float xy = x2*y;
        float xz = x2*z;
        float y2 = y*2;
        float yy = y2*y;
        float yz = y2*z;
        float zz = 2*z*z;
        this.m00 = 1 - yy - zz;
        this.m01 = xy - wz;
        this.m02 = xz + wy;

        this.m10 = xy + wz;
        this.m11 = 1 - xx - zz;
        this.m12 = yz - wx;

        this.m20 = xz - wy;
        this.m21 = yz + wx;
        this.m22 = 1 - xx - yy;
        return this;
    }

    /**
     * r = a * rotate-x(angle)
     * @param a matrix to multiply rotation
     * @param angle angle of rotation
     */
    public Matrix3f rotateX(Matrix3f a, float angle) {
        float s = (float)Math.sin(angle);
        float c = (float)Math.cos(angle);
        float t1, t2;

        t1 = a.m01;
        t2 = a.m02;
        this.m01 = t1 * c + t2 * s;
        this.m02 = t2 * c - t1 * s;

        t1 = a.m11;
        t2 = a.m12;
        this.m11 = t1 * c + t2 * s;
        this.m12 = t2 * c - t1 * s;

        t1 = a.m21;
        t2 = a.m22;
        this.m21 = t1 * c + t2 * s;
        this.m22 = t2 * c - t1 * s;

        this.m00 = a.m00;
        this.m10 = a.m10;
        this.m20 = a.m20;
        return this;
    }

    public Matrix3f rotateY(Matrix3f a, float angle) {
        float s = (float)Math.sin(angle);
        float c = (float)Math.cos(angle);
        float t0 = a.m00;
        float t2 = a.m02;
        this.m00 = t0 * c - t2 * s;
        this.m02 = t0 * s + t2 * c;
        t0 = a.m10;
        t2 = a.m12;
        this.m10 = t0 * c - t2 * s;
        this.m12 = t0 * s + t2 * c;
        t0 = a.m20;
        t2 = a.m22;
        this.m20 = t0 * c - t2 * s;
        this.m22 = t0 * s + t2 * c;
        this.m01 = a.m01;
        this.m11 = a.m11;
        this.m21 = a.m21;
        return this;
    }

    public Matrix3f rotateZ(Matrix3f a, float angle) {
        float s = (float)Math.sin(angle);
        float c = (float)Math.cos(angle);
        float t0 = a.m00;
        float t1 = a.m01;
        this.m00 = t0 * c + t1 * s;
        this.m01 = t1 * c - t0 * s;
        t0 = a.m10;
        t1 = a.m11;
        this.m10 = t0 * c + t1 * s;
        this.m11 = t1 * c - t0 * s;
        t0 = a.m20;
        t1 = a.m21;
        this.m20 = t0 * c + t1 * s;
        this.m21 = t1 * c - t0 * s;
        this.m02 = a.m02;
        this.m12 = a.m12;
        this.m22 = a.m22;
        return this;
    }

    public Matrix3f rotation(float x, float y, float z, float angle) {
        float len = x*x + y*y + z*z;
        if (len < 1e-24f) {
            throw new IllegalArgumentException();
        }
        if (Math.abs(len - 1.0f) > 1e-12f) {
            len = 1.0f/(float)Math.sqrt(len);
            x *= len;
            y *= len;
            z *= len;
        }
        float s = (float)Math.sin(angle);
        float c = (float)Math.cos(angle);
        float t = 1.0f - c;
        float xs = x*s;
        float ys = y*s;
        float zs = z*s;
        float xt = x*t;
        float yt = y*t;
        float zt = z*t;
        this.m00 = xt*x + c;
        this.m01 = xt*y - zs;
        this.m02 = xt*z + ys;

        this.m10 = yt*x + zs;
        this.m11 = yt*y + c;
        this.m12 = yt*z - xs;

        this.m20 = zt*x - ys;
        this.m21 = zt*y + xs;
        this.m22 = zt*z + c;
        return this;
    }

    // computes mul(r, a, rotation_xyzangle)
    // where `a` is a matrix
    // and `rotation` is the rotation matrix produced by
    // a call to rotation(_, x,y,z,angle)
    public Matrix3f rot(Matrix3f a, float x, float y, float z, float angle) {
        float len = x*x + y*y + z*z;
        if (len < 1e-24f) {
            throw new IllegalArgumentException();
        }
        if (Math.abs(len - 1.0f) > 1e-12f) {
            len = 1.0f/(float)Math.sqrt(len);
            x *= len;
            y *= len;
            z *= len;
        }
        float s = (float)Math.sin(angle);
        float c = (float)Math.cos(angle);
        float t = 1 - c;
        float xs = x*s;
        float ys = y*s;
        float zs = z*s;
        float xt = x*t;
        float yt = y*t;
        float zt = z*t;

        float a00 = a.m00;
        float a01 = a.m01;
        float a02 = a.m02;
        float a10 = a.m10;
        float a11 = a.m11;
        float a12 = a.m12;
        float a20 = a.m20;
        float a21 = a.m21;
        float a22 = a.m22;

        float b0 = xt*x + c;
        float b1 = yt*x + zs;
        float b2 = zt*x - ys;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;

        b0 = xt*y - zs;
        b1 = yt*y + c;
        b2 = zt*y + xs;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;

        b0 = xt*z + ys;
        b1 = yt*z - xs;
        b2 = zt*z + c;
        this.m02 = a00*b0 + a01*b1 + a02*b2;
        this.m12 = a10*b0 + a11*b1 + a12*b2;
        this.m22 = a20*b0 + a21*b1 + a22*b2;
        return this;
    }

    public Matrix3f rotate(Matrix3f a, float x, float y, float z, float angle) {
        float len = x*x + y*y + z*z;
        if (len < 1e-24f) return this;
        if (Math.abs(len-1.0f) > 1e-12f) {
            len= 1.0f/(float)Math.sqrt(len);
            x *= len;
            y *= len;
            z *= len;
        }
        float s = (float)Math.sin(angle),
            c = (float)Math.cos(angle),
            t = 1 - c,
            xt = x * t,
            yt = y * t,
            zt = z * t,
            xs = x * s,
            ys = y * s,
            zs = z * s,
            b00 = x * xt + c,
            b01 = y * xt + zs,
            b02 = z * xt - ys,
            b10 = x * yt - zs,
            b11 = y * yt + c,
            b12 = z * yt + xs,
            b20 = x * zt + ys,
            b21 = y * zt - xs,
            b22 = z * zt + c,
            t0 = a.m00,
            t1 = a.m01,
            t2 = a.m02;
        // r00 = cos + x^2 (1 - cos)
        this.m00 = t0 * b00 + t1 * b01 + t2 * b02;
        this.m01 = t0 * b10 + t1 * b11 + t2 * b12;
        this.m02 = t0 * b20 + t1 * b21 + t2 * b22;
        t0 = a.m10;
        t1 = a.m11;
        t2 = a.m12;
        this.m10 = t0 * b00 + t1 * b01 + t2 * b02;
        this.m11 = t0 * b10 + t1 * b11 + t2 * b12;
        this.m12 = t0 * b20 + t1 * b21 + t2 * b22;
        t0 = a.m20;
        t1 = a.m21;
        t2 = a.m22;
        this.m20 = t0 * b00 + t1 * b01 + t2 * b02;
        this.m21 = t0 * b10 + t1 * b11 + t2 * b12;
        this.m22 = t0 * b20 + t1 * b21 + t2 * b22;
        return this;
    }

    public Matrix3f mul(Matrix3f a, Matrix3f b) {
        float a00 = a.m00;
        float a01 = a.m01;
        float a02 = a.m02;
        float a10 = a.m10;
        float a11 = a.m11;
        float a12 = a.m12;
        float a20 = a.m20;
        float a21 = a.m21;
        float a22 = a.m22;

        float b0 = b.m00;
        float b1 = b.m10;
        float b2 = b.m20;
        this.m00 = a00 * b0 + a01 * b1 + a02 * b2;
        this.m10 = a10 * b0 + a11 * b1 + a12 * b2;
        this.m20 = a20 * b0 + a21 * b1 + a22 * b2;

        b0 = b.m01;
        b1 = b.m11;
        b2 = b.m21;
        this.m01 = a00 * b0 + a01 * b1 + a02 * b2;
        this.m11 = a10 * b0 + a11 * b1 + a12 * b2;
        this.m21 = a20 * b0 + a21 * b1 + a22 * b2;

        b0 = b.m02;
        b1 = b.m12;
        b2 = b.m22;
        this.m02 = a00 * b0 + a01 * b1 + a02 * b2;
        this.m12 = a10 * b0 + a11 * b1 + a12 * b2;
        this.m22 = a20 * b0 + a21 * b1 + a22 * b2;
        return this;
    }

    public Matrix3f transposeTimes(AffineTransform3f t1, AffineTransform3f t2) {
        final float a00 = t1.m00;
        final float a01 = t1.m10;
        final float a02 = t1.m20;
        final float a10 = t1.m01;
        final float a11 = t1.m11;
        final float a12 = t1.m21;
        final float a20 = t1.m02;
        final float a21 = t1.m12;
        final float a22 = t1.m22;

        float b0 = t2.m00;
        float b1 = t2.m10;
        float b2 = t2.m20;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;

        b0 = t2.m01;
        b1 = t2.m11;
        b2 = t2.m21;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;

        b0 = t2.m02;
        b1 = t2.m12;
        b2 = t2.m22;
        this.m02 = a00*b0 + a01*b1 + a02*b2;
        this.m12 = a10*b0 + a11*b1 + a12*b2;
        this.m22 = a20*b0 + a21*b1 + a22*b2;

        return this;
    }

    public Matrix3f cwiseApply(Matrix3f m, DoubleFunction<Float> f) {
        this.m00 = (float)f.apply(m.m00);
        this.m10 = (float)f.apply(m.m10);
        this.m20 = (float)f.apply(m.m20);
        this.m01 = (float)f.apply(m.m01);
        this.m11 = (float)f.apply(m.m11);
        this.m21 = (float)f.apply(m.m21);
        this.m02 = (float)f.apply(m.m02);
        this.m12 = (float)f.apply(m.m12);
        this.m22 = (float)f.apply(m.m22);
        return this;
    }

    public Matrix3f abs(Matrix3f m) {
        this.m00 = Math.abs(m.m00);
        this.m10 = Math.abs(m.m10);
        this.m20 = Math.abs(m.m20);
        this.m01 = Math.abs(m.m01);
        this.m11 = Math.abs(m.m11);
        this.m21 = Math.abs(m.m21);
        this.m02 = Math.abs(m.m02);
        this.m12 = Math.abs(m.m12);
        this.m22 = Math.abs(m.m22);
        return this;
    }

    public float colDotX(Vec3f v) {
        return m00*v.x + m10*v.y + m20*v.z;
    }

    public float colDotY(Vec3f v) {
        return m01*v.x + m11*v.y + m21*v.z;
    }

    public float colDotZ(Vec3f v) {
        return m02*v.x + m12*v.y + m22*v.z;
    }

    public float rowDotX(Vec3f v) {
        return m00*v.x + m01*v.y + m02*v.z;
    }

    public float rowDotY(Vec3f v) {
        return m10*v.x + m11*v.y + m12*v.z;
    }

    public float rowDotZ(Vec3f v) {
        return m20*v.x + m21*v.y + m22*v.z;
    }

//    /**
//     * Computes the eigen vector and value of a matrix.
//     * This method violates conventions for this class in
//     * that the output goes to the parameters, instead of this.
//     *
//     * @param vout the eigen vectors (in matrix form)
//     * @param dout the eigen values
//     */
//    public void eigen(Matrix3f vout, Vec3f dout) {
//        Matrix3f R = new Matrix3f(this);
//        final int n = 3;
//        float[] b = new float[3];
//        float[] z = new float[3];
//        float[] d = new float[3];
//
//        vout.setIdentity();
//
//        for (int i=0 ; i<n ; ++i) {
//            b[i] = d[i] = R.getValue(i,i);
//            z[i] = 0;
//        }
//
//        int nrot = 0;
//
//        for (int i=0 ; i<50 ; ++i) {
//            int sm = 0;
//            for (int p=0 ; p<n ; ++p)
//                for (int q=p+1 ; q<n ; ++q)
//                    sm += Math.abs(R.getValue(p, q));
//            if (sm == 0.0f) {
//                dout.set(d[0], d[1], d[2]);
//                return;
//            }
//
//            float thresh  = (i < 3) ? 0.2f * sm / (n*n) : 0.0f;
//
//            for (int p = 0 ; p<n ; ++p) {
//                for (int q = p+i ; q < n ; ++q) {
//                    float g = 100.0f * Math.abs(R.getValue(p, q));
//                    if (i > 3 &&
//                        Math.abs(d[p]) + g == Math.abs(d[p]) &&
//                        Math.abs(d[q]) + g == Math.abs(d[q])) {
//                        R.setValue(p, q, 0.0f);
//                    } else if (Math.abs(R.getValue(p, q)) > thresh) {
//                        float h = d[q] - d[p];
//                        if (Math.abs(h) + g == Math.abs(h))
//                            t = R.get(p, q) / h;
//                        else {
//                            theta = 0.5f * h / R.getValue(p, q);
//
//                        }
//                    }
//                }
//            }
//
//        }
//    }



    public static float determinant(
        float m00, float m01, float m02,
        float m10, float m11, float m12,
        float m20, float m21, float m22)
    {
        return m00*(m11*m22 - m12*m21)
               - m01*(m10*m22 - m12*m20)
               + m02*(m10*m21 - m11*m20);
    }

    public float determinant() {
        return m00*(m11*m22 - m12*m21)
               - m01*(m10*m22 - m12*m20)
               + m02*(m10*m21 - m11*m20);
    }

    public boolean isOrthogonal() {
        // Compute the this x transpose(this)
        // note that Mij == Mji, so we only compute one.

        final float m00 = this.m00;
        final float m01 = this.m01;
        final float m02 = this.m02;
        final float m10 = this.m10;
        final float m11 = this.m11;
        final float m12 = this.m12;
        final float m20 = this.m20;
        final float m21 = this.m21;
        final float m22 = this.m22;

        final float r00 = m00*m00 + m10*m10 + m20*m20;
        final float r11 = m01*m01 + m11*m11 + m21*m21;
        final float r22 = m02*m02 + m12*m12 + m22*m22;
        final float r01 = m00*m01 + m10*m11 + m20*m21;
        final float r02 = m00*m02 + m10*m12 + m20*m22;
        final float r12 = m01*m02 + m11*m12 + m21*m22;

        return Math.abs(1 - r00) < EPSILON &&
               Math.abs(1 - r11) < EPSILON &&
               Math.abs(1 - r22) < EPSILON &&
               Math.abs(r01) < EPSILON &&
               Math.abs(r02) < EPSILON &&
               Math.abs(r12) < EPSILON;
    }

    public boolean isSpecialOrthogonal() {
        return Math.abs(1.0f-determinant()) < EPSILON && isOrthogonal();
    }


    public boolean isSymmetric() {
        return m01 == m10 && m02 == m20 && m12 == m21;
    }

    public float dotRowX(Vec3f v) {
        return m00 * v.x + m01 * v.y + m02 * v.z;
    }

    public float dotRowY(Vec3f v) {
        return m10 * v.x + m11 * v.y + m12 * v.z;
    }

    public float dotRowZ(Vec3f v) {
        return m20 * v.x + m21 * v.y + m22 * v.z;
    }

    public float dotColumnX(Vec3f v) {
        return m00 * v.x + m10 * v.y + m20 * v.z;
    }

    public float dotColumnY(Vec3f v) {
        return m01 * v.x + m11 * v.y + m21 * v.z;
    }

    public float dotColumnZ(Vec3f v) {
        return m02 * v.x + m12 * v.y + m22 * v.z;
    }

    public Matrix3f add(Matrix3f m, float v) {
        m00 = m.m00 + v;
        m01 = m.m01 + v;
        m02 = m.m02 + v;
        m10 = m.m10 + v;
        m11 = m.m11 + v;
        m12 = m.m12 + v;
        m20 = m.m20 + v;
        m21 = m.m21 + v;
        m22 = m.m22 + v;
        return this;
    }

    @Override
    public Matrix3f clone() {
        try {
            return (Matrix3f)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    @Override
    public int hashCode() {
        int h = Float.floatToIntBits(m00);
        h = h * 31 + Float.floatToIntBits(m10);
        h = h * 31 + Float.floatToIntBits(m20);
        h = h * 31 + Float.floatToIntBits(m01);
        h = h * 31 + Float.floatToIntBits(m11);
        h = h * 31 + Float.floatToIntBits(m21);
        h = h * 31 + Float.floatToIntBits(m02);
        h = h * 31 + Float.floatToIntBits(m12);
        h = h * 31 + Float.floatToIntBits(m22);
        return h ^ HASH_BASE;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (null == obj || obj.getClass() != getClass()) return false;

        Matrix3f that = (Matrix3f)obj;
        return MathEx.bitEquals(this.m00, that.m00)
               && MathEx.bitEquals(this.m10, that.m10)
               && MathEx.bitEquals(this.m20, that.m20)
               && MathEx.bitEquals(this.m01, that.m01)
               && MathEx.bitEquals(this.m11, that.m11)
               && MathEx.bitEquals(this.m21, that.m21)
               && MathEx.bitEquals(this.m02, that.m02)
               && MathEx.bitEquals(this.m12, that.m12)
               && MathEx.bitEquals(this.m22, that.m22);
    }

    @Override
    public String toString() {
        return m00 + ", " + m01 + ", " + m02 + "\n" +
               m10 + ", " + m11 + ", " + m12 + "\n" +
               m20 + ", " + m21 + ", " + m22;
    }
}
