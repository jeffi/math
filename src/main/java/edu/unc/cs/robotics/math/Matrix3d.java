package edu.unc.cs.robotics.math;

import java.util.function.DoubleFunction;

/**
 * A 3x3 matrix, specialized for representing rotations in 3D.
 *
 * <p>This class's API generally follows the formula of having {@code this}
 * be the result.  To efficiently use this class, pre-allocate all the
 * result (and operand) matrices.</p>
 */
public final class Matrix3d implements DoubleMatrix, Cloneable {

    private static final int HASH_BASE = Matrix3d.class.getName().hashCode();
    private static final double EPSILON = 1e-9;

    private static final long serialVersionUID = -2720018621259089452L;

    public double m00;
    public double m01;
    public double m02;
    public double m10;
    public double m11;
    public double m12;
    public double m20;
    public double m21;
    public double m22;

    public Matrix3d() {
        m00 = 1.0;
        m11 = 1.0;
        m22 = 1.0;
    }

    public Matrix3d(Matrix3d m) {
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

    public Matrix3d(
        double m00, double m01, double m02,
        double m10, double m11, double m12,
        double m20, double m21, double m22)
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
    public Matrix3d(double[] m) {
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
    public Matrix3d setCoeffs(
        double m00, double m01, double m02,
        double m10, double m11, double m12,
        double m20, double m21, double m22)
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

    public Matrix3d setIdentity() {
        this.m00 = 1.0;
        this.m01 = 0.0;
        this.m02 = 0.0;
        this.m10 = 0.0;
        this.m11 = 1.0;
        this.m12 = 0.0;
        this.m20 = 0.0;
        this.m21 = 0.0;
        this.m22 = 1.0;
        return this;
    }

    public Matrix3d fromTransform(AffineTransform3d t) {
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
    public double getCoeff(int r, int c) {
        if (c < 0 || c > 2)
            throw new IndexOutOfBoundsException();

        double v;

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
    public void setCoeff(int r, int c, double v) {
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

    public Matrix3d transpose(Matrix3d a) {
        // temporary for the case r == a
        double tmp;
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

    public Matrix3d transpose() {
        return transpose(this);
    }

    public Matrix3d setEulerRPY(double x, double y, double z) {
        double cx = Math.cos(x),
            sx = Math.sin(x),
            cy = Math.cos(y),
            sy = Math.sin(y),
            cz = Math.cos(z),
            sz = Math.sin(z),
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

    public Matrix3d fromQuaternion(Quaternion4d q) {
        double w2 = q.w*2;
        double x = q.x;
        double y = q.y;
        double z = q.z;
        double wx = w2*x;
        double wy = w2*y;
        double wz = w2*z;
        double x2 = x*2;
        double xx = x2*x;
        double xy = x2*y;
        double xz = x2*z;
        double y2 = y*2;
        double yy = y2*y;
        double yz = y2*z;
        double zz = 2*z*z;
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
    public Matrix3d rotateX(Matrix3d a, double angle) {
        double s = Math.sin(angle);
        double c = Math.cos(angle);
        double t1, t2;

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

    public Matrix3d rotateY(Matrix3d a, double angle) {
        double s = Math.sin(angle);
        double c = Math.cos(angle);
        double t0 = a.m00;
        double t2 = a.m02;
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

    public Matrix3d rotateZ(Matrix3d a, double angle) {
        double s = Math.sin(angle);
        double c = Math.cos(angle);
        double t0 = a.m00;
        double t1 = a.m01;
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

    public Matrix3d rotation(double x, double y, double z, double angle) {
        double len = x*x + y*y + z*z;
        if (len < 1e-24) {
            throw new IllegalArgumentException();
        }
        if (Math.abs(len - 1.0) > 1e-12) {
            len = 1.0/Math.sqrt(len);
            x *= len;
            y *= len;
            z *= len;
        }
        double s = Math.sin(angle);
        double c = Math.cos(angle);
        double t = 1.0 - c;
        double xs = x*s;
        double ys = y*s;
        double zs = z*s;
        double xt = x*t;
        double yt = y*t;
        double zt = z*t;
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
    public Matrix3d rot(Matrix3d a, double x, double y, double z, double angle) {
        double len = x*x + y*y + z*z;
        if (len < 1e-24) {
            throw new IllegalArgumentException();
        }
        if (Math.abs(len - 1.0) > 1e-12) {
            len = 1.0/Math.sqrt(len);
            x *= len;
            y *= len;
            z *= len;
        }
        double s = Math.sin(angle);
        double c = Math.cos(angle);
        double t = 1 - c;
        double xs = x*s;
        double ys = y*s;
        double zs = z*s;
        double xt = x*t;
        double yt = y*t;
        double zt = z*t;

        double a00 = a.m00;
        double a01 = a.m01;
        double a02 = a.m02;
        double a10 = a.m10;
        double a11 = a.m11;
        double a12 = a.m12;
        double a20 = a.m20;
        double a21 = a.m21;
        double a22 = a.m22;

        double b0 = xt*x + c;
        double b1 = yt*x + zs;
        double b2 = zt*x - ys;
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

    public Matrix3d rotate(Matrix3d a, double x, double y, double z, double angle) {
        double len = x*x + y*y + z*z;
        if (len < 1e-24) return this;
        if (Math.abs(len-1.0) > 1e-12) {
            len= 1.0/Math.sqrt(len);
            x *= len;
            y *= len;
            z *= len;
        }
        double s = Math.sin(angle),
            c = Math.cos(angle),
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

    public Matrix3d mul(Matrix3d a, Matrix3d b) {
        double a00 = a.m00;
        double a01 = a.m01;
        double a02 = a.m02;
        double a10 = a.m10;
        double a11 = a.m11;
        double a12 = a.m12;
        double a20 = a.m20;
        double a21 = a.m21;
        double a22 = a.m22;

        double b0 = b.m00;
        double b1 = b.m10;
        double b2 = b.m20;
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

    public Matrix3d transposeTimes(AffineTransform3d t1, AffineTransform3d t2) {
        final double a00 = t1.m00;
        final double a01 = t1.m10;
        final double a02 = t1.m20;
        final double a10 = t1.m01;
        final double a11 = t1.m11;
        final double a12 = t1.m21;
        final double a20 = t1.m02;
        final double a21 = t1.m12;
        final double a22 = t1.m22;

        double b0 = t2.m00;
        double b1 = t2.m10;
        double b2 = t2.m20;
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

    public Matrix3d cwiseApply(Matrix3d m, DoubleFunction<Double> f) {
        this.m00 = f.apply(m.m00);
        this.m10 = f.apply(m.m10);
        this.m20 = f.apply(m.m20);
        this.m01 = f.apply(m.m01);
        this.m11 = f.apply(m.m11);
        this.m21 = f.apply(m.m21);
        this.m02 = f.apply(m.m02);
        this.m12 = f.apply(m.m12);
        this.m22 = f.apply(m.m22);
        return this;
    }

    public Matrix3d abs(Matrix3d m) {
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

    public double colDotX(Vec3d v) {
        return m00*v.x + m10*v.y + m20*v.z;
    }

    public double colDotY(Vec3d v) {
        return m01*v.x + m11*v.y + m21*v.z;
    }

    public double colDotZ(Vec3d v) {
        return m02*v.x + m12*v.y + m22*v.z;
    }

    public double rowDotX(Vec3d v) {
        return m00*v.x + m01*v.y + m02*v.z;
    }

    public double rowDotY(Vec3d v) {
        return m10*v.x + m11*v.y + m12*v.z;
    }

    public double rowDotZ(Vec3d v) {
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
//    public void eigen(Matrix3d vout, Vec3f dout) {
//        Matrix3d R = new Matrix3d(this);
//        final int n = 3;
//        double[] b = new double[3];
//        double[] z = new double[3];
//        double[] d = new double[3];
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
//            if (sm == 0.0) {
//                dout.set(d[0], d[1], d[2]);
//                return;
//            }
//
//            double thresh  = (i < 3) ? 0.2 * sm / (n*n) : 0.0;
//
//            for (int p = 0 ; p<n ; ++p) {
//                for (int q = p+i ; q < n ; ++q) {
//                    double g = 100.0 * Math.abs(R.getValue(p, q));
//                    if (i > 3 &&
//                        Math.abs(d[p]) + g == Math.abs(d[p]) &&
//                        Math.abs(d[q]) + g == Math.abs(d[q])) {
//                        R.setValue(p, q, 0.0);
//                    } else if (Math.abs(R.getValue(p, q)) > thresh) {
//                        double h = d[q] - d[p];
//                        if (Math.abs(h) + g == Math.abs(h))
//                            t = R.get(p, q) / h;
//                        else {
//                            theta = 0.5 * h / R.getValue(p, q);
//
//                        }
//                    }
//                }
//            }
//
//        }
//    }



    public static double determinant(
        double m00, double m01, double m02,
        double m10, double m11, double m12,
        double m20, double m21, double m22)
    {
        return m00*(m11*m22 - m12*m21)
               - m01*(m10*m22 - m12*m20)
               + m02*(m10*m21 - m11*m20);
    }

    public double determinant() {
        return m00*(m11*m22 - m12*m21)
               - m01*(m10*m22 - m12*m20)
               + m02*(m10*m21 - m11*m20);
    }

    public boolean isOrthogonal() {
        // Compute the this x transpose(this)
        // note that Mij == Mji, so we only compute one.

        final double m00 = this.m00;
        final double m01 = this.m01;
        final double m02 = this.m02;
        final double m10 = this.m10;
        final double m11 = this.m11;
        final double m12 = this.m12;
        final double m20 = this.m20;
        final double m21 = this.m21;
        final double m22 = this.m22;

        final double r00 = m00*m00 + m10*m10 + m20*m20;
        final double r11 = m01*m01 + m11*m11 + m21*m21;
        final double r22 = m02*m02 + m12*m12 + m22*m22;
        final double r01 = m00*m01 + m10*m11 + m20*m21;
        final double r02 = m00*m02 + m10*m12 + m20*m22;
        final double r12 = m01*m02 + m11*m12 + m21*m22;

        return Math.abs(1 - r00) < EPSILON &&
               Math.abs(1 - r11) < EPSILON &&
               Math.abs(1 - r22) < EPSILON &&
               Math.abs(r01) < EPSILON &&
               Math.abs(r02) < EPSILON &&
               Math.abs(r12) < EPSILON;
    }

    public boolean isSpecialOrthogonal() {
        return Math.abs(1.0-determinant()) < EPSILON && isOrthogonal();
    }


    public boolean isSymmetric() {
        return m01 == m10 && m02 == m20 && m12 == m21;
    }

    public double dotRowX(Vec3d v) {
        return m00 * v.x + m01 * v.y + m02 * v.z;
    }

    public double dotRowY(Vec3d v) {
        return m10 * v.x + m11 * v.y + m12 * v.z;
    }

    public double dotRowZ(Vec3d v) {
        return m20 * v.x + m21 * v.y + m22 * v.z;
    }

    public double dotColumnX(Vec3d v) {
        return m00 * v.x + m10 * v.y + m20 * v.z;
    }

    public double dotColumnY(Vec3d v) {
        return m01 * v.x + m11 * v.y + m21 * v.z;
    }

    public double dotColumnZ(Vec3d v) {
        return m02 * v.x + m12 * v.y + m22 * v.z;
    }

    public Matrix3d add(Matrix3d m, double v) {
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
    public Matrix3d clone() {
        try {
            return (Matrix3d)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    @Override
    public int hashCode() {
        long h = Double.doubleToLongBits(m00);
        h = h * 31 + Double.doubleToLongBits(m10);
        h = h * 31 + Double.doubleToLongBits(m20);
        h = h * 31 + Double.doubleToLongBits(m01);
        h = h * 31 + Double.doubleToLongBits(m11);
        h = h * 31 + Double.doubleToLongBits(m21);
        h = h * 31 + Double.doubleToLongBits(m02);
        h = h * 31 + Double.doubleToLongBits(m12);
        h = h * 31 + Double.doubleToLongBits(m22);
        return (int)h ^ (int)(h >>> 32) ^ HASH_BASE;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj) return true;
        if (null == obj || obj.getClass() != getClass()) return false;

        Matrix3d that = (Matrix3d)obj;
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
