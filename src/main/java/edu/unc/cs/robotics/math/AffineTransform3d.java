package edu.unc.cs.robotics.math;

import java.util.Random;

/**
 * An 3D affine transform with doubles.
 * This class represents the transform as a 4x4 transform matrix
 * with the bottom row implied as 0 0 0 1.
 *
 * <pre>
 *   Rxx Rxy Rxz Tx
 *   Ryx Ryy Ryz Ty
 *   Rzx Rzy Rzz Tz
 *   0   0   0   1
 * </pre>
 *
 * This class's API generally follows the formula of having {@code this}
 * be the result.  To efficiently use this class, pre-allocate all the
 * result (and operand) matrices.
 */
public final class AffineTransform3d implements DoubleMatrix, Cloneable {
    private static final int HASH_BASE = AffineTransform3d.class.getName().hashCode();

    private static final long serialVersionUID = -8792057646113184732L;

    public final int ROWS = 4;
    public final int COLUMNS = 4;
    
    public double m00;
    public double m01;
    public double m02;
    public double m03;
    public double m10;
    public double m11;
    public double m12;
    public double m13;
    public double m20;
    public double m21;
    public double m22;
    public double m23;

    /**
     * Creates a new identity transform
     *
     * @see #identity()
     */
    public AffineTransform3d() {
        m00 = 1.0;
        m11 = 1.0;
        m22 = 1.0;
    }

    /**
     * Creates a matrix with the specified coefficients.  No checking
     * on the resulting matrix is performed, thus the produced matrix
     * may not be a valid transformation.  Arguments are provided in
     * row major order (thus source order visually represents the
     * resulting matrix).  In our nomenclature, we use 0-based indexes
     * for the coefficients.
     *
     * @param a00 coeff for row 0, col 0
     * @param a01 coeff for row 0, col 1
     * @param a02 coeff for row 0, col 2
     * @param a03 coeff for row 0, col 3
     * @param a10 coeff for row 1, col 0
     * @param a11 coeff for row 1, col 1
     * @param a12 coeff for row 1, col 2
     * @param a13 coeff for row 1, col 3
     * @param a20 coeff for row 2, col 0
     * @param a21 coeff for row 2, col 1
     * @param a22 coeff for row 2, col 2
     * @param a23 coeff for row 2, col 3
     */
    public AffineTransform3d(
        double a00, double a01, double a02, double a03,
        double a10, double a11, double a12, double a13,
        double a20, double a21, double a22, double a23)
    {
        m00 = a00;
        m01 = a01;
        m02 = a02;
        m03 = a03;
        m10 = a10;
        m11 = a11;
        m12 = a12;
        m13 = a13;
        m20 = a20;
        m21 = a21;
        m22 = a22;
        m23 = a23;
    }

    /**
     * Copy constructor
     *
     * @param src source of copy
     */
    public AffineTransform3d(AffineTransform3d src) {
        set(src);
    }

    /**
     * Creates a transform from a rotation and translation
     *
     * @param q the rotation
     * @param tx the translation in x
     * @param ty the translation in y
     * @param tz the translation in z
     */
    public AffineTransform3d(Quaternion4d q, double tx, double ty, double tz) {
        set(q, tx, ty, tz);
    }

    /**
     * Creates a transform from a rotation and translation
     *
     * @param q the rotation
     * @param t the translation
     */
    public AffineTransform3d(Quaternion4d q, Vec3d t) {
        set(q, t.x, t.y, t.z);
    }

    /**
     * Creates a transform from a rotation and translation
     *
     * @param r the rotation
     * @param x the translation in x
     * @param y the translation in y
     * @param z the translation in z
     */
    public AffineTransform3d(Matrix3d r, double x, double y, double z) {
        set(r, x, y, z);
    }

    /**
     * Creates a transform from a rotation and translation
     *
     * @param r the rotation
     * @param t the translation
     */
    public AffineTransform3d(Matrix3d r, Vec3d t) {
        set(r, t.x, t.y, t.z);
    }

    // ================================================================
    // Set methods (methods that change the entire matrix)
    // ================================================================

    /**
     * Sets this matrix to the identity matrix.
     *
     * @return {@code this}
     */
    public AffineTransform3d identity() {
        m00 = 1.0;
        m01 = 0.0;
        m02 = 0.0;
        m03 = 0.0;
        m10 = 0.0;
        m11 = 1.0;
        m12 = 0.0;
        m13 = 0.0;
        m20 = 0.0;
        m21 = 0.0;
        m22 = 1.0;
        m23 = 0.0;
        return this;
    }

    /**
     * Sets this transform to be a copy of the argument.
     *
     * @param src the matrix to copy
     * @return {@code this}
     */
    public AffineTransform3d set(AffineTransform3d src) {
        m00 = src.m00;
        m01 = src.m01;
        m02 = src.m02;
        m03 = src.m03;
        m10 = src.m10;
        m11 = src.m11;
        m12 = src.m12;
        m13 = src.m13;
        m20 = src.m20;
        m21 = src.m21;
        m22 = src.m22;
        m23 = src.m23;
        return this;
    }

    /**
     * Sets the rotation and translation to match the arguments.
     *
     * @param q the rotation
     * @param tx the translation in x
     * @param ty the translation in y
     * @param tz the translation in z
     * @return {@code this}
     */
    public AffineTransform3d set(Quaternion4d q, double tx, double ty, double tz) {
        final double w2 = q.w*2.0;
        final double x = q.x;
        final double y = q.y;
        final double z = q.z;
        final double wx = w2*x;
        final double wy = w2*y;
        final double wz = w2*z;
        final double x2 = x*2.0;
        final double xx = x2*x;
        final double xy = x2*y;
        final double xz = x2*z;
        final double y2 = y*2.0;
        final double yy = y2*y;
        final double yz = y2*z;
        final double zz = 2.0*z*z;

        m00 = 1.0 - yy - zz;
        m01 = xy - wz;
        m02 = xz + wy;

        m10 = xy + wz;
        m11 = 1.0 - xx - zz;
        m12 = yz - wx;

        m20 = xz - wy;
        m21 = yz + wx;
        m22 = 1.0 - xx - yy;

        m03 = tx;
        m13 = ty;
        m23 = tz;
        return this;
    }


    /**
     * Sets the rotation and translation to match the arguments.
     *
     * @param q the rotation to set
     * @param t the translation to set
     * @return this
     */
    public AffineTransform3d set(Quaternion4d q, Vec3d t) {
        return set(q, t.x, t.y, t.z);
    }

    public AffineTransform3d set(Matrix3d r, double x, double y, double z) {
        m00 = r.m00;
        m01 = r.m01;
        m02 = r.m02;
        m10 = r.m10;
        m11 = r.m11;
        m12 = r.m12;
        m20 = r.m20;
        m21 = r.m21;
        m22 = r.m22;
        m03 = x;
        m13 = y;
        m23 = z;
        return this;
    }

    public AffineTransform3d set(Matrix3d r, Vec3d t) {
        return set(r, t.x, t.y, t.z);
    }

    /**
     * Sets this matrix to be a rotation matrix around an arbitrary axis.
     * The axis does not need to be normalized.  If the axis is zero, the result
     * will be a matrix with a lot of NaNs.
     *
     * @param x the x value of the axis
     * @param y the y value of the axis
     * @param z the z value of the axis
     * @param angle the angle of rotation (in radians)
     * @return {@code this}
     */
    public AffineTransform3d rotation(double x, double y, double z, double angle) {
        double len = x*x + y*y + z*z;
        if (Math.abs(1.0 - len) > 1e-9) {
            len = Math.sqrt(len);
            x/=len;
            y/=len;
            z/=len;
        }

        final double c = Math.cos(angle);
        final double s = Math.sin(angle);
        m03 = 0.0;
        m13 = 0.0;
        m23 = 0.0;

        final double t = 1.0 - c;
        final double zs = z*s;
        final double ys = y*s;
        final double xs = x*s;
        final double xt = x*t;
        m00 = x*xt + c;
        m10 = y*xt + zs;
        m20 = z*xt - ys;

        final double yt = y*t;
        m01 = x*yt - zs;
        m11 = y*yt + c;
        m21 = z*yt + xs;

        final double zt = z*t;
        m02 = x*zt + ys;
        m12 = y*zt - xs;
        m22 = z*zt + c;

        return this;
    }

    public AffineTransform3d rotation(Vec3d axis, double angle) {
        return rotation(axis.x, axis.y, axis.z, angle);
    }

    /**
     * Sets this transform to be a rotation about the x-axis.
     *
     * @param angle the angle (in radians)
     * @return {@code this}
     */
    public AffineTransform3d rotationX(double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);
        m03 = 0.0;
        m13 = 0.0;
        m23 = 0.0;

        m00 = 1.0;
        m10 = 0.0;
        m20 = 0.0;

        m01 = 0.0;
        m11 = c;
        m21 = s;

        m02 = 0.0;
        m12 = -s;
        m22 = c;

        return this;
    }

    /**
     * Sets this transform to be a rotation about the y-axis.
     *
     * @param angle the angle (in radians)
     * @return {@code this}
     */
    public AffineTransform3d rotationY(double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);
        m03 = 0.0;
        m13 = 0.0;
        m23 = 0.0;

        m00 = c;
        m10 = 0.0;
        m20 = -s;

        m01 = 0.0;
        m11 = 1.0;
        m21 = 0.0;

        m02 = s;
        m12 = 0.0;
        m22 = c;

        return this;
    }

    /**
     * Sets this transform to be a rotation about the z-axis.
     *
     * @param angle the angle (in radians)
     * @return {@code this}
     */
    public AffineTransform3d rotationZ(double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);
        m03 = 0.0;
        m13 = 0.0;
        m23 = 0.0;

        m00 = c;
        m10 = s;
        m20 = 0.0;

        m01 = -s;
        m11 = c;
        m21 = 0.0;

        m02 = 0.0;
        m12 = 0.0;
        m22 = 1.0;

        return this;
    }

    /**
     * Sets this transform to be a rotation transform.  If setting the translation
     * as well, use {@link #set(Quaternion4d, Vec3d)}.
     *
     * @param q the rotation to use
     * @return {@code this}
     * @see #set(Quaternion4d, Vec3d)
     */
    public AffineTransform3d rotation(Quaternion4d q) {
        return set(q, 0.0, 0.0, 0.0);
    }

    /**
     * Sets a the transform to represent an Euler roll-pitch-yaw (XYZ) transform.
     * This is equivalent to the matrix created by:
     *
     * <pre>
     * t.setIdentity();
     * t.rotateZ(yaw);
     * t.rotateY(pitch);
     * t.rotateX(roll);
     * </pre>
     *
     * @param roll rotation about X
     * @param pitch rotation about Y
     * @param yaw rotation about Z
     * @return this
     */
    public AffineTransform3d rotationRPY(double roll, double pitch, double yaw) {
        final double cz = Math.cos(yaw);
        final double sz = Math.sin(yaw);
        final double cy = Math.cos(pitch);
        final double sy = Math.sin(pitch);
        final double cx = Math.cos(roll);
        final double sx = Math.sin(roll);
        final double czsy = cz*sy;
        final double szsy = sz*sy;

        m03 = 0.0;
        m13 = 0.0;
        m23 = 0.0;

        m00 = cz*cy;
        m10 = sz*cy;
        m20 = -sy;

        m01 = czsy*sx - sz*cx;
        m11 = cz*cx + szsy*sx;
        m21 = cy*sx;

        m02 = czsy*cx + sz*sx;
        m12 = szsy*cx - cz*sx;
        m22 = cy*cx;

        return this;
    }

    /**
     * Sets this transform to be a translation in the form:
     *
     * <pre>
     *     1 0 0 x
     *     0 1 0 y
     *     0 0 1 z
     * </pre>
     *
     * @param x the x translation
     * @param y the y translation
     * @param z the z translation
     * @return {@code this}
     */
    public AffineTransform3d translation(double x, double y, double z) {
        m00 = 1.0;
        m01 = 0.0;
        m02 = 0.0;
        m03 = x;
        m10 = 0.0;
        m11 = 1.0;
        m12 = 0.0;
        m13 = y;
        m20 = 0.0;
        m21 = 0.0;
        m22 = 1.0;
        m23 = z;
        return this;
    }

    public AffineTransform3d translation(Vec3d v) {
        return translation(v.x, v.y, v.z);
    }

    /**
     * Equivalent to {@code setTranslate(x, y, z).mulRPY(roll, pitch, yaw)}
     *
     * @param x translation x
     * @param y translation y
     * @param z translation z
     * @param roll roll angle
     * @param pitch pitch angle
     * @param yaw yaw angle
     * @return this
     */
    public AffineTransform3d setTranslationRPY(
        double x, double y, double z,
        double roll, double pitch, double yaw)
    {
        final double cz = Math.cos(yaw);
        final double sz = Math.sin(yaw);
        final double cy = Math.cos(pitch);
        final double sy = Math.sin(pitch);
        final double cx = Math.cos(roll);
        final double sx = Math.sin(roll);
        final double czsy = cz*sy;
        final double szsy = sz*sy;

        m03 = x;
        m13 = y;
        m23 = z;

        m00 = cz*cy;
        m10 = sz*cy;
        m20 = -sy;

        m01 = czsy*sx - sz*cx;
        m11 = cz*cx + szsy*sx;
        m21 = cy*sx;

        m02 = czsy*cx + sz*sx;
        m12 = szsy*cx - cz*sx;
        m22 = cy*cx;

        return this;
    }

    // ================================================================
    // Interface Methods
    // ================================================================

    @Override
    public int rows() {
        return ROWS;
    }

    @Override
    public int columns() {
        return COLUMNS;
    }

    @Override
    public int size() {
        return ROWS*COLUMNS;
    }

    @Override
    public double getCoeff(int r, int c) {
        if (c < 0 || c >= COLUMNS)
            throw new MatrixIndexOutOfBoundsException(this, r, c);

        switch (r*COLUMNS + c) {
        case 0: return m00;
        case 1: return m01;
        case 2: return m02;
        case 3: return m03;
        case 4: return m10;
        case 5: return m11;
        case 6: return m12;
        case 7: return m13;
        case 8: return m20;
        case 9: return m21;
        case 10: return m22;
        case 11: return m23;
        case 12: return 0.0;
        case 13: return 0.0;
        case 14: return 0.0;
        case 15: return 1.0;
        default: throw new MatrixIndexOutOfBoundsException(this, r, c);
        }
    }

    @Override
    public void setCoeff(int r, int c, double value) {
        if (c < 0 || c >= COLUMNS)
            throw new MatrixIndexOutOfBoundsException(this, r, c);

        switch (r*COLUMNS + c) {
        case 0: m00 = value; break;
        case 1: m01 = value; break;
        case 2: m02 = value; break;
        case 3: m03 = value; break;
        case 4: m10 = value; break;
        case 5: m11 = value; break;
        case 6: m12 = value; break;
        case 7: m13 = value; break;
        case 8: m20 = value; break;
        case 9: m21 = value; break;
        case 10: m22 = value; break;
        case 11: m23 = value; break;
        default: throw new MatrixIndexOutOfBoundsException(this, r, c);
        }
    }

    // ================================================================
    // Multiplication methods
    // These method compute A x B, where A and/or B are varying
    // transformation forms
    // ================================================================

    /**
     * Computes {@code this = a * b}.
     *
     * @param a the first operand
     * @param b the second operand
     * @return {@code this}
     */
    public AffineTransform3d mul(AffineTransform3d a, AffineTransform3d b) {
        final double a00 = a.m00;
        final double a01 = a.m01;
        final double a02 = a.m02;
        final double a10 = a.m10;
        final double a11 = a.m11;
        final double a12 = a.m12;
        final double a20 = a.m20;
        final double a21 = a.m21;
        final double a22 = a.m22;
        double b0 = b.m00;
        double b1 = b.m10;
        double b2 = b.m20;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;
        b0 = b.m01;
        b1 = b.m11;
        b2 = b.m21;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;
        b0 = b.m02;
        b1 = b.m12;
        b2 = b.m22;
        this.m02 = a00*b0 + a01*b1 + a02*b2;
        this.m12 = a10*b0 + a11*b1 + a12*b2;
        this.m22 = a20*b0 + a21*b1 + a22*b2;
        b0 = b.m03;
        b1 = b.m13;
        b2 = b.m23;
        this.m03 = a00*b0 + a01*b1 + a02*b2 + a.m03;
        this.m13 = a10*b0 + a11*b1 + a12*b2 + a.m13;
        this.m23 = a20*b0 + a21*b1 + a22*b2 + a.m23;
        return this;
    }

    /**
     * Multiplies this matrix by the argument.  @{code this x b}
     *
     * @param b the right-hand side operand of the multiplication.
     * @return {@code this}
     */
    public AffineTransform3d mul(AffineTransform3d b) {
        return mul(this, b);
    }


    /**
     * Rotates the argument transform about the specified axis by the specified angle.
     * The axis does not need to be normalized prior to calling.
     * Stores the result in {@code this}
     *
     * @param a the transform to rotate
     * @param x the x component of the axis
     * @param y the y component of the axis
     * @param z the z component of the axis
     * @param angle the rotation in radians
     * @throws IllegalArgumentException if the axis is (0,0,0)
     * @return {@code this}
     */
    public AffineTransform3d rotate(AffineTransform3d a, double x, double y, double z, double angle) {
        double len = x*x + y*y + z*z;
        if (Math.abs(1.0 - len) > 1e-9) {
            len = Math.sqrt(len);
            if (len == 0.0) {
                throw new IllegalArgumentException("0 axis");
            }
            x /= len;
            y /= len;
            z /= len;
        }

        final double c = Math.cos(angle);
        final double s = Math.sin(angle);

        final double t = 1.0 - c;
        final double zs = z*s;
        final double ys = y*s;
        final double xs = x*s;

        final double a00 = a.m00;
        final double a01 = a.m01;
        final double a02 = a.m02;
        final double a10 = a.m10;
        final double a11 = a.m11;
        final double a12 = a.m12;
        final double a20 = a.m20;
        final double a21 = a.m21;
        final double a22 = a.m22;

        final double xt = x*t;
        double b0 = x*xt + c;
        double b1 = y*xt + zs;
        double b2 = z*xt - ys;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;
        final double yt = y*t;
        b0 = x*yt - zs;
        b1 = y*yt + c;
        b2 = z*yt + xs;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;
        final double zt = z*t;
        b0 = x*zt + ys;
        b1 = y*zt - xs;
        b2 = z*zt + c;
        this.m02 = a00*b0 + a01*b1 + a02*b2;
        this.m12 = a10*b0 + a11*b1 + a12*b2;
        this.m22 = a20*b0 + a21*b1 + a22*b2;
        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;
        return this;
    }

    public AffineTransform3d rotate(AffineTransform3d a, Vec3d axis, double angle) {
        return rotate(a, axis.x, axis.y, axis.z, angle);
    }

    /**
     * Rotates this matrix by computing {@code this = this * rotation(x,y,z,angle)}.
     *
     * @param x the x value of the axis
     * @param y the y value of the axis
     * @param z the z value of the axis
     * @param angle the angle of rotation (in radians)
     * @return {@code this}
     */
    public AffineTransform3d rotate(double x, double y, double z, double angle) {
        return rotate(this, x, y, z, angle);
    }

    public AffineTransform3d rotate(Vec3d axis, double angle) {
        return rotate(this, axis.x, axis.y, axis.z, angle);
    }


    /**
     * Computes the product of {@code a} and a rotation matrix
     * around the x-axis.
     *
     * @param a left-hand operand of multiplication
     * @param angle the angle of rotation about the x-axis
     * @return {@code this}
     */
    public AffineTransform3d rotateX(AffineTransform3d a, double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);

        final double a01 = a.m01;
        final double a11 = a.m11;
        final double a21 = a.m21;
        final double a02 = a.m02;
        final double a12 = a.m12;
        final double a22 = a.m22;

        this.m00 = a.m00;
        this.m10 = a.m10;
        this.m20 = a.m20;

        this.m01 = a01*c + a02*s;
        this.m11 = a11*c + a12*s;
        this.m21 = a21*c + a22*s;

        this.m02 = a02*c - a01*s;
        this.m12 = a12*c - a11*s;
        this.m22 = a22*c - a21*s;

        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;

        return this;
    }

    /**
     * Rotates this matrix about the x-axis by computing
     * {@code this = this * rotationX(angle)}
     *
     * @param angle the angle (in radians)
     * @return {@code this}
     */
    public AffineTransform3d rotateX(double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);

        final double a01 = m01;
        final double a11 = m11;
        final double a21 = m21;
        final double a02 = m02;
        final double a12 = m12;
        final double a22 = m22;

        this.m01 = a01*c + a02*s;
        this.m11 = a11*c + a12*s;
        this.m21 = a21*c + a22*s;

        this.m02 = a02*c - a01*s;
        this.m12 = a12*c - a11*s;
        this.m22 = a22*c - a21*s;

        return this;
    }


    public AffineTransform3d rotateY(AffineTransform3d a, double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);

        final double a00 = a.m00;
        final double a10 = a.m10;
        final double a20 = a.m20;
        final double a02 = a.m02;
        final double a12 = a.m12;
        final double a22 = a.m22;

        this.m00 = a00*c - a02*s;
        this.m10 = a10*c - a12*s;
        this.m20 = a20*c - a22*s;

        this.m01 = a.m01;
        this.m11 = a.m11;
        this.m21 = a.m21;

        this.m02 = a00*s + a02*c;
        this.m12 = a10*s + a12*c;
        this.m22 = a20*s + a22*c;

        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;

        return this;
    }

    public AffineTransform3d rotateY(double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);

        final double a00 = m00;
        final double a10 = m10;
        final double a20 = m20;
        final double a02 = m02;
        final double a12 = m12;
        final double a22 = m22;

        this.m00 = a00*c - a02*s;
        this.m10 = a10*c - a12*s;
        this.m20 = a20*c - a22*s;

        this.m02 = a00*s + a02*c;
        this.m12 = a10*s + a12*c;
        this.m22 = a20*s + a22*c;

        return this;
    }

    public AffineTransform3d rotateZ(AffineTransform3d a, double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);

        final double a00 = a.m00;
        final double a10 = a.m10;
        final double a20 = a.m20;
        final double a01 = a.m01;
        final double a11 = a.m11;
        final double a21 = a.m21;

        this.m00 = a00*c + a01*s;
        this.m10 = a10*c + a11*s;
        this.m20 = a20*c + a21*s;

        this.m01 = a01*c - a00*s;
        this.m11 = a11*c - a10*s;
        this.m21 = a21*c - a20*s;

        this.m02 = a.m02;
        this.m12 = a.m12;
        this.m22 = a.m22;

        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;

        return this;
    }

    public AffineTransform3d rotateZ(double angle) {
        final double c = Math.cos(angle);
        final double s = Math.sin(angle);

        final double a00 = m00;
        final double a10 = m10;
        final double a20 = m20;
        final double a01 = m01;
        final double a11 = m11;
        final double a21 = m21;

        this.m00 = a00*c + a01*s;
        this.m10 = a10*c + a11*s;
        this.m20 = a20*c + a21*s;

        this.m01 = a01*c - a00*s;
        this.m11 = a11*c - a10*s;
        this.m21 = a21*c - a20*s;

        return this;
    }

    /**
     * Rotates a transform by a rotation matrix.
     *
     * @param a the transform to rotate
     * @param r the rotation to apply
     * @return {@code this}
     */
    public AffineTransform3d rotate(AffineTransform3d a, Matrix3d r) {
        final double a00 = a.m00;
        final double a01 = a.m01;
        final double a02 = a.m02;
        final double a10 = a.m10;
        final double a11 = a.m11;
        final double a12 = a.m12;
        final double a20 = a.m20;
        final double a21 = a.m21;
        final double a22 = a.m22;
        double b0 = r.m00;
        double b1 = r.m10;
        double b2 = r.m20;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;
        b0 = r.m01;
        b1 = r.m11;
        b2 = r.m21;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;
        b0 = r.m02;
        b1 = r.m12;
        b2 = r.m22;
        this.m02 = a00*b0 + a01*b1 + a02*b2;
        this.m12 = a10*b0 + a11*b1 + a12*b2;
        this.m22 = a20*b0 + a21*b1 + a22*b2;
        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;
        return this;
    }

    /**
     * Equivalent to {@code a.rotate(a, r)}
     * @param r the rotation
     * @return {@code this}
     * @see #rotate(AffineTransform3d, Matrix3d)
     */
    public AffineTransform3d rotate(Matrix3d r) {
        return rotate(this, r);
    }

    /**
     * Rotates a transform by a quaternion and stores the result in {@code this}.
     *
     * @param a the transform to rotate
     * @param q the rotation
     * @return {@code this}
     */
    public AffineTransform3d rotate(AffineTransform3d a, Quaternion4d q) {
        final double w2 = q.w*2.0;
        final double x = q.x;
        final double y = q.y;
        final double z = q.z;
        final double wx = w2*x;
        final double wy = w2*y;
        final double wz = w2*z;
        final double x2 = x*2.0;
        final double xx = x2*x;
        final double xy = x2*y;
        final double xz = x2*z;
        final double y2 = y*2.0;
        final double yy = y2*y;
        final double yz = y2*z;
        final double zz = 2.0*z*z;

        final double a00 = a.m00;
        final double a01 = a.m01;
        final double a02 = a.m02;
        final double a10 = a.m10;
        final double a11 = a.m11;
        final double a12 = a.m12;
        final double a20 = a.m20;
        final double a21 = a.m21;
        final double a22 = a.m22;
        double b0 = 1.0 - yy - zz;
        double b1 = xy + wz;
        double b2 = xz - wy;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;
        b0 = xy - wz;
        b1 = 1.0 - xx - zz;
        b2 = yz + wx;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;
        b0 = xz + wy;
        b1 = yz - wx;
        b2 = 1.0 - xx - yy;
        this.m02 = a00*b0 + a01*b1 + a02*b2;
        this.m12 = a10*b0 + a11*b1 + a12*b2;
        this.m22 = a20*b0 + a21*b1 + a22*b2;
        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;
        return this;
    }

    public AffineTransform3d rotate(Quaternion4d q) {
        return rotate(this, q);
    }

    public AffineTransform3d rotateRPY(AffineTransform3d a, double r, double p, double y) {
        final double cz = Math.cos(y);
        final double sz = Math.sin(y);
        final double cy = Math.cos(p);
        final double sy = Math.sin(p);
        final double cx = Math.cos(r);
        final double sx = Math.sin(r);
        final double czsy = cz*sy;
        final double szsy = sz*sy;

        final double a00 = a.m00;
        final double a01 = a.m01;
        final double a02 = a.m02;
        final double a10 = a.m10;
        final double a11 = a.m11;
        final double a12 = a.m12;
        final double a20 = a.m20;
        final double a21 = a.m21;
        final double a22 = a.m22;
        double b0 = cz*cy;
        double b1 = sz*cy;
        double b2 = -sy;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;
        b0 = czsy*sx - sz*cx;
        b1 = cz*cx + szsy*sx;
        b2 = cy*sx;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;
        b0 = czsy*cx + sz*sx;
        b1 = szsy*cx - cz*sx;
        b2 = cy*cx;
        this.m02 = a00*b0 + a01*b1 + a02*b2;
        this.m12 = a10*b0 + a11*b1 + a12*b2;
        this.m22 = a20*b0 + a21*b1 + a22*b2;
        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;
        return this;
    }

    public AffineTransform3d rotateRPY(double r, double p, double y) {
        return rotateRPY(this, r, p, y);
    }

    public AffineTransform3d setScale(double x, double y, double z) {
        m10 = m20 = m01 = m21 = m02 = m12 = m03 = m13 = m23 = 0.0;
        m00 = x;
        m11 = y;
        m22 = z;
        return this;
    }

    public AffineTransform3d scale(AffineTransform3d a, double x, double y, double z) {
        this.m00 = a.m00*x;
        this.m10 = a.m10*x;
        this.m20 = a.m20*x;

        this.m01 = a.m01*y;
        this.m11 = a.m11*y;
        this.m21 = a.m21*y;

        this.m02 = a.m02*z;
        this.m12 = a.m12*z;
        this.m22 = a.m22*z;

        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;

        return this;
    }

    public AffineTransform3d scale(double x, double y, double z) {
        this.m00 *= x;
        this.m10 *= x;
        this.m20 *= x;

        this.m01 *= y;
        this.m11 *= y;
        this.m21 *= y;

        this.m02 *= z;
        this.m12 *= z;
        this.m22 *= z;

        return this;
    }


    /**
     * Computes the multiplication of a transform by a translation and stores
     * the result in {@code this}.
     *
     * <pre>
     *                | 1 0 0 x |
     *     this = a * | 0 1 0 y |
     *                | 0 0 1 z |
     * </pre>
     * @param a the matrix to multiply
     * @param x the translation in x
     * @param y the translation in y
     * @param z the translation in z
     * @return {@code this}
     */
    public AffineTransform3d translate(AffineTransform3d a, double x, double y, double z) {
        final double a00 = a.m00;
        final double a10 = a.m10;
        final double a20 = a.m20;
        final double a01 = a.m01;
        final double a11 = a.m11;
        final double a21 = a.m21;
        final double a02 = a.m02;
        final double a12 = a.m12;
        final double a22 = a.m22;

        this.m00 = a00;
        this.m10 = a10;
        this.m20 = a20;

        this.m01 = a01;
        this.m11 = a11;
        this.m21 = a21;

        this.m02 = a02;
        this.m12 = a12;
        this.m22 = a22;

        this.m03 = a00*x + a01*y + a02*z + a.m03;
        this.m13 = a10*x + a11*y + a12*z + a.m13;
        this.m23 = a20*x + a21*y + a22*z + a.m23;

        return this;
    }

    public AffineTransform3d translate(AffineTransform3d a, Vec3d t) {
        return translate(a, t.x, t.y, t.z);
    }

    /**
     * Multiplies this matrix by a translation matrix, and returns this.
     * Equivalent to {@code mul(createTranslate(x,y,z))}, however since the
     * structure of the translation matrix is known, the number of operations
     * is greatly reduced.
     *
     * <pre>
     *                | 1 0 0 x |
     *  this = this x | 0 1 0 y |
     *                | 0 0 1 z |
     * </pre>
     *
     * @param x the translation in x
     * @param y the translation in y
     * @param z the translation in z
     * @return {@code this}
     */
    public AffineTransform3d translate(double x, double y, double z) {
        // inlined translate(this, x, y, z) call:
        this.m03 += m00*x + m01*y + m02*z;
        this.m13 += m10*x + m11*y + m12*z;
        this.m23 += m20*x + m21*y + m22*z;
        return this;
    }

    public AffineTransform3d translate(Vec3d t) {
        return translate(t.x, t.y, t.z);
    }


    /**
     * Sets the rotation to be a random rotation.
     * Leaves the translation untouched.
     *
     * @param rand source of randomness
     * @return this
     */
    public AffineTransform3d randomRotation(Random rand) {
        // The following code is equivalent to
        // Quaterniond q = new Quaterniond();
        // q.random(rand);
        // setRotation(q);
        // setTranslation( ... );
        // return this;

        // 2 trig pairs
        // 2 sqrts
        // 19 mul
        // 13 add

//        double x0 = rand.nextDouble();
//        double r1 = Math.sqrt(1.0 - x0);
//        double r2 = Math.sqrt(x0);
//        double t1 = rand.nextDouble()*(Math.PI*2);
//        double t2 = rand.nextDouble()*(Math.PI*2);
//        double c1 = Math.cos(t1);
//        double s1 = Math.sin(t1);
//        double c2 = Math.cos(t2);
//        double s2 = Math.sin(t2);
//        double qw = s1 * r1;
//        double qx = c1 * r1;
//        double qy = s2 * r2;
//        double qz = c2 * r2;
//        double w2 = qw*2;
//        double wx = w2*qx;
//        double wy = w2*qy;
//        double wz = w2*qz;
//        double x2 = qx*2;
//        double xx = x2*qx;
//        double xy = x2*qy;
//        double xz = x2*qz;
//        double y2 = qy*2;
//        double yy = y2*qy;
//        double yz = y2*qz;
//        double zz = 2*qz*qz;
//        m00 = 1 - yy - zz;
//        m01 = xy - wz;
//        m02 = xz + wy;
//
//        m10 = xy + wz;
//        m11 = 1 - xx - zz;
//        m12 = yz - wx;
//
//        m20 = xz - wy;
//        m21 = yz + wx;
//        m22 = 1 - xx - yy;




        // 2 trig pairs
        // 1 sqrt
        // 16 mul
        // 12 add
        double r2r2 = rand.nextDouble()*2;
        double r1r1 = (2.0 - r2r2);
        double r1r2 = Math.sqrt(r1r1*r2r2);
        double t1 = rand.nextDouble()*(Math.PI*2);
        double t2 = rand.nextDouble()*(Math.PI*2);
        double c1 = Math.cos(t1);
        double s1 = Math.sin(t1);
        double c2 = Math.cos(t2);
        double s2 = Math.sin(t2);
        double c1s2 = c1*s2;
        double c2s1 = c2*s1;
        double c1c2 = c1*c2;
        double s1s2 = s1*s2;
        double r1r1c1c1_1 = 1 - r1r1*c1*c1;
        double r2r2c2s2 = r2r2*c2*s2;
        double r1r1c2s1 = r1r1*c1*s1;
        double r2r2s2s2 = r2r2*s2*s2;

        m00 = 1 - r2r2;
        m11 = r1r1c1c1_1 + r2r2s2s2 - r2r2; // c2*c2;
        m22 = r1r1c1c1_1 - r2r2s2s2;

        m01 = r1r2*(c1s2 - c2s1);
        m10 = r1r2*(c1s2 + c2s1);

        m02 = r1r2*(c1c2 + s1s2);
        m20 = r1r2*(c1c2 - s1s2);

        m12 = r2r2c2s2 - r1r1c2s1;
        m21 = r2r2c2s2 + r1r1c2s1;

        return this;
    }

    /**
     * {@code invert(this)}.
     *
     * @see #inverse(AffineTransform3d)
     * @return this
     */
    public AffineTransform3d inverse() {
        return inverse(this);
    }

    /**
     * Performs a (slow) inversion of this transform.  This inversion handles
     * any invertible matrix (rotation, translation, scale, sheer).  If the
     * matrix is only a rotation and a translation, the rotation can be
     * transposed, and the translation negated.
     *
     * @param m the matrix to invert
     * @return this
     * @throws ArithmeticException if the matrix cannot be inverted
     */
    public AffineTransform3d inverse(AffineTransform3d m) {
        // See http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
        // This code follows that algorithm, though is significantly reduced based
        // upon the bottom row begin fixed at [0,0,0,1]

        double a11 = m.m11;
        double a21 = m.m21;
        double a22 = m.m22;
        double a12 = m.m12;
        double inv00 = a11*a22 - a12*a21;
        double a01 = m.m01;
        double a02 = m.m02;
        double inv01 = a02*a21 - a01*a22;
        double inv02 = a01*a12 - a02*a11;
        double a00 = m.m00;
        double a10 = m.m10;
        double a20 = m.m20;
        double det = a00*inv00 + a10*inv01 + a20*inv02;

        if (det == 0.0)
            throw new ArithmeticException("transform not invertible");

        det = 1.0/det;

        this.m00 = inv00*det;
        this.m01 = inv01*det;
        this.m02 = inv02*det;
        this.m10 = (a12*a20 - a10*a22)*det;
        this.m11 = (a00*a22 - a02*a20)*det;
        this.m12 = (a02*a10 - a00*a12)*det;
        this.m20 = (a10*a21 - a11*a20)*det;
        this.m21 = (a01*a20 - a00*a21)*det;
        this.m22 = (a00*a11 - a01*a10)*det;

        double tx = -m.m03;
        double ty = m.m13;
        double tz = m.m23;
        this.m03 = tx*m00 - ty*m01 - tz*m02;
        this.m13 = tx*m10 - ty*m11 - tz*m12;
        this.m23 = tx*m20 - ty*m21 - tz*m22;

        return this;
    }

    /**
     * Computes the determinant of the affine transform.  Since the bottom
     * row of the implied 4x4 matrix is [0, 0, 0, 1], the determinant is
     * the same as the determinant of just the rotation component of the matrix.
     * If the rotation remains special orthogonal the determinant should always be
     * 1.0, however the determinant is not a sufficient check for special orthogonality.
     *
     * @return the determinant of the transform matrix
     */
    public double determinant() {
        return m00*(m11*m22 - m12*m21)
               + m10*(m02*m21 - m01*m22)
               + m20*(m01*m12 - m02*m11);
    }

    /**
     * Checks if the rotation component of this transform is orthogonal.
     * A matrix is orthogonal if its transpose is its inverse.
     *
     * @return true if the rotation is orthogonal
     */
    public boolean isRotationOrthogonal() {
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

        return Math.abs(1 - r00) < 1e-9 &&
               Math.abs(1 - r11) < 1e-9 &&
               Math.abs(1 - r22) < 1e-9 &&
               Math.abs(r01) < 1e-9 &&
               Math.abs(r02) < 1e-9 &&
               Math.abs(r12) < 1e-9;
    }

    /**
     * Checks if the rotation component of this transform is special orthogonal.
     * A special orthogonal matrix is a proper rotation (no sheers, scales, etc...),
     * and the rotation can be inverted by a transpose.  A matrix is special
     * orthogonal if its is orthogonal and its determinant is 1 (the only other
     * possibility for an orthogonal matrix is -1).
     *
     * <p>In general, the rotation component of an AffineTransform3d instance will
     * remain special orthogonal as long as only rotations (and translations) are
     * applied to the transform.  However numeric precision may cause the matrix
     * to drift out of orthogonality and require fixing.</p>
     *
     * <p>When an AffineTransform3d instance has a special orthogonal rotation
     * component, one can use the fastRigidInverse methods.</p>
     *
     * @return true if the rotation is special orthogonal
     */
    public boolean isRotationSpecialOrthogonal() {
        return Math.abs(1.0 - determinant()) < 1e-9 && isRotationOrthogonal();
    }


    /**
     * Computes an inverse transform with the assumption that
     * the argument has only rotations and translations.
     *
     * @param m the transform to invert
     * @return {@code this}
     */
    public AffineTransform3d fastRigidInverse(AffineTransform3d m) {
        double tx = -m.m03;
        double ty = m.m13;
        double tz = m.m23;
        // in case this == m, store these in temporaries
        double tmp10 = m.m01;
        double tmp20 = m.m02;
        double tmp21 = m.m12;
        m03 = tx*(m00 = m.m00) - ty*(m01 = m.m10) - tz*(m02 = m.m20);
        m13 = tx*(m10 = tmp10) - ty*(m11 = m.m11) - tz*(m12 = m.m21);
        m23 = tx*(m20 = tmp20) - ty*(m21 = tmp21) - tz*(m22 = m.m22);
        return this;
    }

    /**
     * Equivalent to:
     * <pre>
     *     x.fastRigidInverse(x);
     * </pre>
     *
     * @see #fastRigidInverse(AffineTransform3d)
     * @return {@code this}
     */
    public AffineTransform3d fastRigidInverse() {
        return fastRigidInverse(this);
    }

    /**
     * Sets the translation to a random value in the provided
     * range.  Leaves the rotation untouched.
     *
     * @param rand source of randomness
     * @param xMin minimum value for x
     * @param xMax maximum value for x
     * @param yMin minimum value for y
     * @param yMax maximum value for y
     * @param zMin minimum value for z
     * @param zMax maximum value for z
     * @return this
     */
    public AffineTransform3d randomTranslation(
        Random rand,
        double xMin, double xMax,
        double yMin, double yMax,
        double zMin, double zMax)
    {
        this.m03 = rand.nextDouble() * (xMax - xMin) + xMin;
        this.m13 = rand.nextDouble() * (yMax - yMin) + yMin;
        this.m23 = rand.nextDouble() * (zMax - zMin) + zMin;
        return this;
    }

    /**
     * Sets the rotation and translation to random values.
     *
     * @param rand source of randomness
     * @param xMin minimum value for x
     * @param xMax maximum value for x
     * @param yMin minimum value for y
     * @param yMax maximum value for y
     * @param zMin minimum value for z
     * @param zMax maximum value for z
     * @return this
     */
    public AffineTransform3d random(
        Random rand,
        double xMin, double xMax,
        double yMin, double yMax,
        double zMin, double zMax)
    {
        return randomRotation(rand)
            .randomTranslation(rand, xMin, xMax, yMin, yMax, zMin, zMax);
    }

    /**
     * Computes t1'*t2 considering only the rotation.  The result is stored
     * in this.  The result's translation component is set to 0.
     *
     * @param t1 transposed operand
     * @param t2 operand
     * @return this
     */
    public AffineTransform3d rotationTransposeTimes(AffineTransform3d t1, AffineTransform3d t2) {
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
        this.m00 = a00 * b0 + a01 * b1 + a02 * b2;
        this.m10 = a10 * b0 + a11 * b1 + a12 * b2;
        this.m20 = a20 * b0 + a21 * b1 + a22 * b2;

        b0 = t2.m01;
        b1 = t2.m11;
        b2 = t2.m21;
        this.m01 = a00 * b0 + a01 * b1 + a02 * b2;
        this.m11 = a10 * b0 + a11 * b1 + a12 * b2;
        this.m21 = a20 * b0 + a21 * b1 + a22 * b2;

        b0 = t2.m02;
        b1 = t2.m12;
        b2 = t2.m22;
        this.m02 = a00 * b0 + a01 * b1 + a02 * b2;
        this.m12 = a10 * b0 + a11 * b1 + a12 * b2;
        this.m22 = a20 * b0 + a21 * b1 + a22 * b2;

        this.m03 = 0;
        this.m13 = 0;
        this.m23 = 0;

        return this;
    }

    /**
     * Generates a coordinate system given a single direction vector
     *
     * @param w a vector
     * @return this
     */
    public AffineTransform3d generateCoordinateSystem(Vec3d w) {
        m00 = w.x;
        m10 = w.y;
        m20 = w.z;

        if (Math.abs(w.x) >= Math.abs(w.y)) {
            double invLen = 1.0/Math.sqrt(m00*m00 + m20*m20);
            m01 = -m20*invLen;
            m11 = 0;
            m21 = m00*invLen;
            m02 = m10*m21;
            m12 = m20*m01 - m00*m21;
            m22 = -m10*m01;
        } else {
            double inv_length = 1.0/Math.sqrt(m10*m10 + m20*m20);
            m01 = 0;
            m11 = m20*inv_length;
            m21 = -m10*inv_length;
            m02 = m10*m21 - m20*m11;
            m12 = -m00*m21;
            m22 = m00*m11;
        }

        m03 = 0;
        m13 = 0;
        m23 = 0;

        return this;
    }

    public AffineTransform3d setAxes(Vec3d u, Vec3d v, Vec3d w) {
        m00 = u.x;
        m10 = u.y;
        m20 = u.z;

        m01 = v.x;
        m11 = v.y;
        m21 = v.z;

        m02 = w.x;
        m12 = w.y;
        m22 = w.z;

        return this;
    }

    public AffineTransform3d makeRotationOrthogonal() {
        return makeRotationOrthogonal(this);
    }

    public AffineTransform3d makeRotationOrthogonal(AffineTransform3d t) {
        double ux = t.m00;
        double uy = t.m10;
        double uz = t.m20;
        double vx = t.m01;
        double vy = t.m11;
        double vz = t.m21;

        // compute w = u x v
        double wx = uy * vz - uz * vy;
        double wy = uz * vx - ux * vz;
        double wz = ux * vy - uy * vx;

        // normalize w
        double dw = Math.sqrt(wx*wx + wy*wy + wz*wz);
        wx /= dw;
        wy /= dw;
        wz /= dw;

        // normalize v
        double dv = Math.sqrt(vx*vx + vy*vy + vz*vz);
        vx /= dv;
        vy /= dv;
        vz /= dv;

        // compute u = v x w
        m00 = vy * wz - vz * wy;
        m10 = vz * wx - vx * wz;
        m20 = vx * wy - vy * wx;

        m01 = vx;
        m11 = vy;
        m21 = vz;

        m02 = wx;
        m12 = wy;
        m22 = wz;

        return this;
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

    public double dotColumn(int axis, Vec3d v) {
        switch (axis) {
        case 0: return dotColumnX(v);
        case 1: return dotColumnY(v);
        case 2: return dotColumnZ(v);
        default:
            throw new IndexOutOfBoundsException();
        }

    }

    /**
     * Compares two transforms for mathematical equality within an
     * error bound.  If precise equality is desired, while treating
     * {@code +0.0f == -0.0f} and {@code NaN != NaN}, then use this
     * method with {@code epsilon = 0.0f}.  (Contrast with {@link #equals(Object)},
     * which treats {@code +0.0f != -0.0f} and {@code NaN == NaN}).
     *
     * @param that the other transform to compare
     * @param epsilon the error range (must be >= 0)
     * @return true if each element is within epsilon of the other.
     */
    public boolean equals(AffineTransform3d that, double epsilon) {
        if (epsilon < 0) {
            throw new IllegalArgumentException();
        }

        return Math.abs(this.m00 - that.m00) <= epsilon
               && Math.abs(this.m01 - that.m01) <= epsilon
               && Math.abs(this.m02 - that.m02) <= epsilon
               && Math.abs(this.m03 - that.m03) <= epsilon
               && Math.abs(this.m10 - that.m10) <= epsilon
               && Math.abs(this.m11 - that.m11) <= epsilon
               && Math.abs(this.m12 - that.m12) <= epsilon
               && Math.abs(this.m13 - that.m13) <= epsilon
               && Math.abs(this.m20 - that.m20) <= epsilon
               && Math.abs(this.m21 - that.m21) <= epsilon
               && Math.abs(this.m22 - that.m22) <= epsilon
               && Math.abs(this.m23 - that.m23) <= epsilon;
    }

    @Override
    public AffineTransform3d clone() {
        try {
            return (AffineTransform3d)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    /**
     * Produces a hashcode for this transform.  Compatible with {@link #equals(Object)}.
     *
     * @return the hash code
     */
    @Override
    public int hashCode() {
        long h =   Double.doubleToLongBits(m00);
        h = h*31 + Double.doubleToLongBits(m01);
        h = h*31 + Double.doubleToLongBits(m02);
        h = h*31 + Double.doubleToLongBits(m03);
        h = h*31 + Double.doubleToLongBits(m10);
        h = h*31 + Double.doubleToLongBits(m11);
        h = h*31 + Double.doubleToLongBits(m12);
        h = h*31 + Double.doubleToLongBits(m13);
        h = h*31 + Double.doubleToLongBits(m20);
        h = h*31 + Double.doubleToLongBits(m21);
        h = h*31 + Double.doubleToLongBits(m22);
        h = h*31 + Double.doubleToLongBits(m23);
        return (int)h ^ (int)(h >>> 32) ^ HASH_BASE;
    }

    /**
     * Equals method compatible with hashing.  Avoid using this method for comparing
     * two transforms, as it is likely not what you want.  In order to be compatible
     * with hashing contract, the coefficients are compared by their bit values, this
     * means that unlike '==', that NaN's are considered equivalent, and 0.0 is not
     * the same as -0.0.
     *
     * @param o the object to compare
     * @return true if objects are bit-wise equals.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || o.getClass() != getClass()) return false;
        final AffineTransform3d that = (AffineTransform3d)o;
        return Double.doubleToLongBits(this.m00) == Double.doubleToLongBits(that.m00)
            && Double.doubleToLongBits(this.m01) == Double.doubleToLongBits(that.m01)
            && Double.doubleToLongBits(this.m02) == Double.doubleToLongBits(that.m02)
            && Double.doubleToLongBits(this.m03) == Double.doubleToLongBits(that.m03)
            && Double.doubleToLongBits(this.m10) == Double.doubleToLongBits(that.m10)
            && Double.doubleToLongBits(this.m11) == Double.doubleToLongBits(that.m11)
            && Double.doubleToLongBits(this.m12) == Double.doubleToLongBits(that.m12)
            && Double.doubleToLongBits(this.m13) == Double.doubleToLongBits(that.m13)
            && Double.doubleToLongBits(this.m20) == Double.doubleToLongBits(that.m20)
            && Double.doubleToLongBits(this.m21) == Double.doubleToLongBits(that.m21)
            && Double.doubleToLongBits(this.m22) == Double.doubleToLongBits(that.m22)
            && Double.doubleToLongBits(this.m23) == Double.doubleToLongBits(that.m23);
    }

    /**
     * Returns a string with 1 line per row of the transform matrix (not including the
     * implied [0, 0, 0 1] row.
     *
     * @return a string with 1 line per row of the transform matrix.
     */
    @Override
    public String toString() {
        return m00 + ", " + m01 + ", " + m02 + ", " + m03 + "\n"
             + m10 + ", " + m11 + ", " + m12 + ", " + m13 + "\n"
             + m20 + ", " + m21 + ", " + m22 + ", " + m23;
    }
}
