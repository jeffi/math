// Automatically generated from AffineTransform3d.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import java.util.Random;

/**
 * An 3D affine transform with floats.
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
public final class AffineTransform3f implements FloatMatrix, Cloneable {
    private static final int HASH_BASE = AffineTransform3f.class.getName().hashCode();

    private static final long serialVersionUID = -8792057646113184732L;

    public final int ROWS = 4;
    public final int COLUMNS = 4;
    
    public float m00;
    public float m01;
    public float m02;
    public float m03;
    public float m10;
    public float m11;
    public float m12;
    public float m13;
    public float m20;
    public float m21;
    public float m22;
    public float m23;

    /**
     * Creates a new identity transform
     *
     * @see #identity()
     */
    public AffineTransform3f() {
        m00 = 1.0f;
        m11 = 1.0f;
        m22 = 1.0f;
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
    public AffineTransform3f(
        float a00, float a01, float a02, float a03,
        float a10, float a11, float a12, float a13,
        float a20, float a21, float a22, float a23)
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
    public AffineTransform3f(AffineTransform3f src) {
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
    public AffineTransform3f(Quaternion4f q, float tx, float ty, float tz) {
        set(q, tx, ty, tz);
    }

    /**
     * Creates a transform from a rotation and translation
     *
     * @param q the rotation
     * @param t the translation
     */
    public AffineTransform3f(Quaternion4f q, Vec3f t) {
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
    public AffineTransform3f(Matrix3f r, float x, float y, float z) {
        set(r, x, y, z);
    }

    /**
     * Creates a transform from a rotation and translation
     *
     * @param r the rotation
     * @param t the translation
     */
    public AffineTransform3f(Matrix3f r, Vec3f t) {
        set(r, t.x, t.y, t.z);
    }

    /**
     * Creates a new transform that is a translation with the specified
     * values for translation.
     *
     * @param tx translation in x
     * @param ty translation in y
     * @param tz translation in z
     */
    public AffineTransform3f(float tx, float ty, float tz) {
        setTranslation(tx, ty, tz);
    }

    /**
     * Creates a new transform that is a translation with the specified
     * vector for its translation.
     *
     * @param t the translation in (x,y,z)
     */
    public AffineTransform3f(Vec3f t) {
        setTranslation(t.x, t.y, t.z);
    }

    /**
     * Creates a new transform that with the values from an array with
     * the 12 elements in column-major order.
     *
     * @param m the values to set for this transform matrix.
     * @see #set(float[])
     */
    public AffineTransform3f(float[] m) {
        set(m);
    }

    /**
     * Creates a new transform that with the values from an array with
     * the 12 elements in column-major order.
     *
     * @param m the values to set for this transform matrix.
     * @param offset the offset of the (0, 0) coefficient
     * @param stride the offset between columns
     * @see #set(float[], int, int)
     */
    public AffineTransform3f(float[] m, int offset, int stride) {
        set(m, offset, stride);
    }

    // ================================================================
    // Set methods (methods that change the entire matrix)
    // ================================================================

    /**
     * Sets this matrix to the identity matrix.
     *
     * @return {@code this}
     */
    public AffineTransform3f identity() {
        m00 = 1.0f;
        m01 = 0.0f;
        m02 = 0.0f;
        m03 = 0.0f;
        m10 = 0.0f;
        m11 = 1.0f;
        m12 = 0.0f;
        m13 = 0.0f;
        m20 = 0.0f;
        m21 = 0.0f;
        m22 = 1.0f;
        m23 = 0.0f;
        return this;
    }

    /**
     * Sets this transform to be a copy of the argument.
     *
     * @param src the matrix to copy
     * @return {@code this}
     */
    public AffineTransform3f set(AffineTransform3f src) {
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
    public AffineTransform3f set(Quaternion4f q, float tx, float ty, float tz) {
        final float w2 = q.w*2.0f;
        final float x = q.x;
        final float y = q.y;
        final float z = q.z;
        final float wx = w2*x;
        final float wy = w2*y;
        final float wz = w2*z;
        final float x2 = x*2.0f;
        final float xx = x2*x;
        final float xy = x2*y;
        final float xz = x2*z;
        final float y2 = y*2.0f;
        final float yy = y2*y;
        final float yz = y2*z;
        final float zz = 2.0f*z*z;

        m00 = 1.0f - yy - zz;
        m01 = xy - wz;
        m02 = xz + wy;

        m10 = xy + wz;
        m11 = 1.0f - xx - zz;
        m12 = yz - wx;

        m20 = xz - wy;
        m21 = yz + wx;
        m22 = 1.0f - xx - yy;

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
     * @return {@code this}
     */
    public AffineTransform3f set(Quaternion4f q, Vec3f t) {
        return set(q, t.x, t.y, t.z);
    }

    /**
     * Sets this transform to be a rotation and translation.
     *
     * @param r the rotation
     * @param tx translation in x
     * @param ty translation in y
     * @param tz translation in z
     * @return {@code this}
     */
    public AffineTransform3f set(Matrix3f r, float tx, float ty, float tz) {
        m00 = r.m00;
        m01 = r.m01;
        m02 = r.m02;
        m10 = r.m10;
        m11 = r.m11;
        m12 = r.m12;
        m20 = r.m20;
        m21 = r.m21;
        m22 = r.m22;
        m03 = tx;
        m13 = ty;
        m23 = tz;
        return this;
    }

    /**
     * Sets this transform to be a rotation and translation.
     *
     * @param r the rotation
     * @param t the translation
     * @return {@code this}
     */
    public AffineTransform3f set(Matrix3f r, Vec3f t) {
        return set(r, t.x, t.y, t.z);
    }

    /**
     * Sets the values of this from a 3x4 matrix stored as a column-major array.
     *
     * @param m the matrix to set
     * @return {@code this}
     * @see #AffineTransform3f(float[])
     */
    public AffineTransform3f set(float[] m) {
        m00 = m[0];
        m10 = m[1];
        m20 = m[2];
        m01 = m[3];
        m11 = m[4];
        m21 = m[5];
        m02 = m[6];
        m12 = m[7];
        m22 = m[8];
        m03 = m[9];
        m13 = m[10];
        m23 = m[11];
        return this;
    }

    /**
     * Sets the value of this from a 3x4 matrix stored as a column-major array.
     * The values are pulled from the specified offset and stride within the
     * array.
     *
     * @param m the matrix to set
     * @param offset the offset of the (0, 0) element
     * @param stride the offset between columns
     * @return {@code this}
     * @see #AffineTransform3f(float[], int, int)
     */
    public AffineTransform3f set(float[] m, int offset, int stride) {
        m00 = m[offset];
        m10 = m[offset + 1];
        m20 = m[offset + 2];
        offset += stride;
        m01 = m[offset];
        m11 = m[offset + 1];
        m21 = m[offset + 2];
        offset += stride;
        m02 = m[offset];
        m12 = m[offset + 1];
        m22 = m[offset + 2];
        offset += stride;
        m03 = m[offset];
        m13 = m[offset + 1];
        m23 = m[offset + 2];
        return this;
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
    public AffineTransform3f rotation(float x, float y, float z, float angle) {
        float len = x*x + y*y + z*z;
        if (Math.abs(1.0f - len) > 1e-9f) {
            len = (float)Math.sqrt(len);
            x/=len;
            y/=len;
            z/=len;
        }

        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);
        m03 = 0.0f;
        m13 = 0.0f;
        m23 = 0.0f;

        final float t = 1.0f - c;
        final float zs = z*s;
        final float ys = y*s;
        final float xs = x*s;
        final float xt = x*t;
        m00 = x*xt + c;
        m10 = y*xt + zs;
        m20 = z*xt - ys;

        final float yt = y*t;
        m01 = x*yt - zs;
        m11 = y*yt + c;
        m21 = z*yt + xs;

        final float zt = z*t;
        m02 = x*zt + ys;
        m12 = y*zt - xs;
        m22 = z*zt + c;

        return this;
    }

    public AffineTransform3f rotation(Vec3f axis, float angle) {
        return rotation(axis.x, axis.y, axis.z, angle);
    }

    /**
     * Sets this transform to be a rotation about the x-axis.
     *
     * @param angle the angle (in radians)
     * @return {@code this}
     */
    public AffineTransform3f rotationX(float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);
        m03 = 0.0f;
        m13 = 0.0f;
        m23 = 0.0f;

        m00 = 1.0f;
        m10 = 0.0f;
        m20 = 0.0f;

        m01 = 0.0f;
        m11 = c;
        m21 = s;

        m02 = 0.0f;
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
    public AffineTransform3f rotationY(float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);
        m03 = 0.0f;
        m13 = 0.0f;
        m23 = 0.0f;

        m00 = c;
        m10 = 0.0f;
        m20 = -s;

        m01 = 0.0f;
        m11 = 1.0f;
        m21 = 0.0f;

        m02 = s;
        m12 = 0.0f;
        m22 = c;

        return this;
    }

    /**
     * Sets this transform to be a rotation about the z-axis.
     *
     * @param angle the angle (in radians)
     * @return {@code this}
     */
    public AffineTransform3f rotationZ(float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);
        m03 = 0.0f;
        m13 = 0.0f;
        m23 = 0.0f;

        m00 = c;
        m10 = s;
        m20 = 0.0f;

        m01 = -s;
        m11 = c;
        m21 = 0.0f;

        m02 = 0.0f;
        m12 = 0.0f;
        m22 = 1.0f;

        return this;
    }

    /**
     * Sets this transform to be a rotation transform.  If setting the translation
     * as well, use {@link #set(Quaternion4f, Vec3f)}.
     *
     * @param q the rotation to use
     * @return {@code this}
     * @see #set(Quaternion4f, Vec3f)
     */
    public AffineTransform3f rotation(Quaternion4f q) {
        return set(q, 0.0f, 0.0f, 0.0f);
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
    public AffineTransform3f rotationRPY(float roll, float pitch, float yaw) {
        final float cz = (float)Math.cos(yaw);
        final float sz = (float)Math.sin(yaw);
        final float cy = (float)Math.cos(pitch);
        final float sy = (float)Math.sin(pitch);
        final float cx = (float)Math.cos(roll);
        final float sx = (float)Math.sin(roll);
        final float czsy = cz*sy;
        final float szsy = sz*sy;

        m03 = 0.0f;
        m13 = 0.0f;
        m23 = 0.0f;

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
    public AffineTransform3f translation(float x, float y, float z) {
        m00 = 1.0f;
        m01 = 0.0f;
        m02 = 0.0f;
        m03 = x;
        m10 = 0.0f;
        m11 = 1.0f;
        m12 = 0.0f;
        m13 = y;
        m20 = 0.0f;
        m21 = 0.0f;
        m22 = 1.0f;
        m23 = z;
        return this;
    }

    public AffineTransform3f translation(Vec3f v) {
        return translation(v.x, v.y, v.z);
    }

    /**
     * Equivalent to {@code translation(x, y, z).rotateRPY(roll, pitch, yaw)}
     *
     * @param x translation x
     * @param y translation y
     * @param z translation z
     * @param roll roll angle
     * @param pitch pitch angle
     * @param yaw yaw angle
     * @return this
     */
    public AffineTransform3f translationRPY(
        float x, float y, float z,
        float roll, float pitch, float yaw)
    {
        final float cz = (float)Math.cos(yaw);
        final float sz = (float)Math.sin(yaw);
        final float cy = (float)Math.cos(pitch);
        final float sy = (float)Math.sin(pitch);
        final float cx = (float)Math.cos(roll);
        final float sx = (float)Math.sin(roll);
        final float czsy = cz*sy;
        final float szsy = sz*sy;

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
    // Set translation independent of rotation
    // ================================================================

    /**
     * Sets the translation to the specified values while leaving the
     * rotation untouched.
     *
     * @param x the new value for the translation in x
     * @param y the new value for the translation in y
     * @param z the new value for the translation in z
     * @return {@code this}
     */
    public AffineTransform3f setTranslation(float x, float y, float z) {
        m03 = x;
        m13 = y;
        m23 = z;
        return this;
    }

    /**
     * Sets the translation to the specified values while leaving the
     * rotation untouched.
     *
     * @param t the new value for the translation
     * @return {@code this}
     * @see #setTranslation(Vec3f)
     */
    public AffineTransform3f setTranslation(Vec3f t) {
        return setTranslation(t.x, t.y, t.z);
    }

    /**
     * Adds the specified value to the translation.  This is equivalent
     * to multiplying this transform by a translation, with the translation
     * as the left-hand side operand.
     *
     * <pre>
     *     a.preTranslate(x, y, z);
     *     // same as:
     *     // a.mul(new AffineTransform3f().translation(x, y, z), a);
     * </pre>
     *
     * @param x the translation in x to add
     * @param y the translation in y to add
     * @param z the translation in z to add
     * @return this
     */
    public AffineTransform3f preTranslate(float x, float y, float z) {
        this.m03 += x;
        this.m13 += y;
        this.m23 += z;
        return this;
    }

    /**
     * Adds the specified vector to the translation.
     *
     * @param t the translation to add
     * @return this
     * @see #preTranslate(float, float, float)
     */
    public AffineTransform3f preTranslate(Vec3f t) {
        return preTranslate(t.x, t.y, t.z);
    }


    // ================================================================
    // Interface Methods
    // ================================================================

    @Override
    public int rows() {
        return 4;
    }

    /**
     * Returns the number of columns in the underlying matrix.
     * In this case the return value is 4 since it includes the
     * implied last row.
     *
     * @return 4
     */
    @Override
    public int columns() {
        return 4;
    }

    @Override
    public int size() {
        return 16;
    }

    @Override
    public float getCoeff(int r, int c) {
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
        case 12: return 0.0f;
        case 13: return 0.0f;
        case 14: return 0.0f;
        case 15: return 1.0f;
        default: throw new MatrixIndexOutOfBoundsException(this, r, c);
        }
    }

    @Override
    public void setCoeff(int r, int c, float value) {
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
        case 12:
        case 13:
        case 14:
        case 15:
            throw new UnsupportedOperationException(
                "last row is not modifiable");
        default:
            throw new MatrixIndexOutOfBoundsException(this, r, c);
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
    public AffineTransform3f mul(AffineTransform3f a, AffineTransform3f b) {
        final float a00 = a.m00;
        final float a01 = a.m01;
        final float a02 = a.m02;
        final float a10 = a.m10;
        final float a11 = a.m11;
        final float a12 = a.m12;
        final float a20 = a.m20;
        final float a21 = a.m21;
        final float a22 = a.m22;
        float b0 = b.m00;
        float b1 = b.m10;
        float b2 = b.m20;
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

    // Alternate multiplication formulation that is slightly slower
//    public AffineTransform3f mul2(AffineTransform3f a, AffineTransform3f b) {
//        final float b00 = b.m00;
//        final float b10 = b.m10;
//        final float b20 = b.m20;
//
//        final float b01 = b.m01;
//        final float b11 = b.m11;
//        final float b21 = b.m21;
//
//        final float b02 = b.m02;
//        final float b12 = b.m12;
//        final float b22 = b.m22;
//
//        final float b03 = b.m03;
//        final float b13 = b.m13;
//        final float b23 = b.m23;
//
//        float a0 = a.m00;
//        float a1 = a.m01;
//        float a2 = a.m02;
//        this.m00 = a0*b00 + a1*b10 + a2*b20;
//        this.m01 = a0*b01 + a1*b11 + a2*b21;
//        this.m02 = a0*b02 + a1*b12 + a2*b22;
//        this.m03 = a0*b03 + a1*b13 + a2*b23 + a.m03;
//
//        a0 = a.m10;
//        a1 = a.m11;
//        a2 = a.m12;
//        this.m10 = a0*b00 + a1*b10 + a2*b20;
//        this.m11 = a0*b01 + a1*b11 + a2*b21;
//        this.m12 = a0*b02 + a1*b12 + a2*b22;
//        this.m13 = a0*b03 + a1*b13 + a2*b23 + a.m13;
//
//        a0 = a.m20;
//        a1 = a.m21;
//        a2 = a.m22;
//        this.m20 = a0*b00 + a1*b10 + a2*b20;
//        this.m21 = a0*b01 + a1*b11 + a2*b21;
//        this.m22 = a0*b02 + a1*b12 + a2*b22;
//        this.m23 = a0*b03 + a1*b13 + a2*b23 + a.m23;
//        return this;
//    }

    /**
     * Multiplies this matrix by the argument.  @{code this x b}
     *
     * @param b the right-hand side operand of the multiplication.
     * @return {@code this}
     */
    public AffineTransform3f mul(AffineTransform3f b) {
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
    public AffineTransform3f rotate(AffineTransform3f a, float x, float y, float z, float angle) {
        float len = x*x + y*y + z*z;
        if (Math.abs(1.0f - len) > 1e-9f) {
            len = (float)Math.sqrt(len);
            if (len == 0.0f) {
                throw new IllegalArgumentException("0 axis");
            }
            x /= len;
            y /= len;
            z /= len;
        }

        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);

        final float t = 1.0f - c;
        final float zs = z*s;
        final float ys = y*s;
        final float xs = x*s;

        final float a00 = a.m00;
        final float a01 = a.m01;
        final float a02 = a.m02;
        final float a10 = a.m10;
        final float a11 = a.m11;
        final float a12 = a.m12;
        final float a20 = a.m20;
        final float a21 = a.m21;
        final float a22 = a.m22;

        final float xt = x*t;
        float b0 = x*xt + c;
        float b1 = y*xt + zs;
        float b2 = z*xt - ys;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;
        final float yt = y*t;
        b0 = x*yt - zs;
        b1 = y*yt + c;
        b2 = z*yt + xs;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;
        final float zt = z*t;
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

    public AffineTransform3f rotate(AffineTransform3f a, Vec3f axis, float angle) {
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
    public AffineTransform3f rotate(float x, float y, float z, float angle) {
        return rotate(this, x, y, z, angle);
    }

    public AffineTransform3f rotate(Vec3f axis, float angle) {
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
    public AffineTransform3f rotateX(AffineTransform3f a, float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);

        final float a01 = a.m01;
        final float a11 = a.m11;
        final float a21 = a.m21;
        final float a02 = a.m02;
        final float a12 = a.m12;
        final float a22 = a.m22;

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
    public AffineTransform3f rotateX(float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);

        final float a01 = m01;
        final float a11 = m11;
        final float a21 = m21;
        final float a02 = m02;
        final float a12 = m12;
        final float a22 = m22;

        this.m01 = a01*c + a02*s;
        this.m11 = a11*c + a12*s;
        this.m21 = a21*c + a22*s;

        this.m02 = a02*c - a01*s;
        this.m12 = a12*c - a11*s;
        this.m22 = a22*c - a21*s;

        return this;
    }


    public AffineTransform3f rotateY(AffineTransform3f a, float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);

        final float a00 = a.m00;
        final float a10 = a.m10;
        final float a20 = a.m20;
        final float a02 = a.m02;
        final float a12 = a.m12;
        final float a22 = a.m22;

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

    public AffineTransform3f rotateY(float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);

        final float a00 = m00;
        final float a10 = m10;
        final float a20 = m20;
        final float a02 = m02;
        final float a12 = m12;
        final float a22 = m22;

        this.m00 = a00*c - a02*s;
        this.m10 = a10*c - a12*s;
        this.m20 = a20*c - a22*s;

        this.m02 = a00*s + a02*c;
        this.m12 = a10*s + a12*c;
        this.m22 = a20*s + a22*c;

        return this;
    }

    public AffineTransform3f rotateZ(AffineTransform3f a, float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);

        final float a00 = a.m00;
        final float a10 = a.m10;
        final float a20 = a.m20;
        final float a01 = a.m01;
        final float a11 = a.m11;
        final float a21 = a.m21;

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

    public AffineTransform3f rotateZ(float angle) {
        final float c = (float)Math.cos(angle);
        final float s = (float)Math.sin(angle);

        final float a00 = m00;
        final float a10 = m10;
        final float a20 = m20;
        final float a01 = m01;
        final float a11 = m11;
        final float a21 = m21;

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
    public AffineTransform3f rotate(AffineTransform3f a, Matrix3f r) {
        final float a00 = a.m00;
        final float a01 = a.m01;
        final float a02 = a.m02;
        final float a10 = a.m10;
        final float a11 = a.m11;
        final float a12 = a.m12;
        final float a20 = a.m20;
        final float a21 = a.m21;
        final float a22 = a.m22;
        float b0 = r.m00;
        float b1 = r.m10;
        float b2 = r.m20;
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
     * @see #rotate(AffineTransform3f, Matrix3f)
     */
    public AffineTransform3f rotate(Matrix3f r) {
        return rotate(this, r);
    }

    /**
     * Rotates a transform by a quaternion and stores the result in {@code this}.
     *
     * @param a the transform to rotate
     * @param q the rotation
     * @return {@code this}
     */
    public AffineTransform3f rotate(AffineTransform3f a, Quaternion4f q) {
        final float w2 = q.w*2.0f;
        final float x = q.x;
        final float y = q.y;
        final float z = q.z;
        final float wx = w2*x;
        final float wy = w2*y;
        final float wz = w2*z;
        final float x2 = x*2.0f;
        final float xx = x2*x;
        final float xy = x2*y;
        final float xz = x2*z;
        final float y2 = y*2.0f;
        final float yy = y2*y;
        final float yz = y2*z;
        final float zz = 2.0f*z*z;

        final float a00 = a.m00;
        final float a01 = a.m01;
        final float a02 = a.m02;
        final float a10 = a.m10;
        final float a11 = a.m11;
        final float a12 = a.m12;
        final float a20 = a.m20;
        final float a21 = a.m21;
        final float a22 = a.m22;
        float b0 = 1.0f - yy - zz;
        float b1 = xy + wz;
        float b2 = xz - wy;
        this.m00 = a00*b0 + a01*b1 + a02*b2;
        this.m10 = a10*b0 + a11*b1 + a12*b2;
        this.m20 = a20*b0 + a21*b1 + a22*b2;
        b0 = xy - wz;
        b1 = 1.0f - xx - zz;
        b2 = yz + wx;
        this.m01 = a00*b0 + a01*b1 + a02*b2;
        this.m11 = a10*b0 + a11*b1 + a12*b2;
        this.m21 = a20*b0 + a21*b1 + a22*b2;
        b0 = xz + wy;
        b1 = yz - wx;
        b2 = 1.0f - xx - yy;
        this.m02 = a00*b0 + a01*b1 + a02*b2;
        this.m12 = a10*b0 + a11*b1 + a12*b2;
        this.m22 = a20*b0 + a21*b1 + a22*b2;
        this.m03 = a.m03;
        this.m13 = a.m13;
        this.m23 = a.m23;
        return this;
    }

    public AffineTransform3f rotate(Quaternion4f q) {
        return rotate(this, q);
    }

    public AffineTransform3f rotateRPY(AffineTransform3f a, float r, float p, float y) {
        final float cz = (float)Math.cos(y);
        final float sz = (float)Math.sin(y);
        final float cy = (float)Math.cos(p);
        final float sy = (float)Math.sin(p);
        final float cx = (float)Math.cos(r);
        final float sx = (float)Math.sin(r);
        final float czsy = cz*sy;
        final float szsy = sz*sy;

        final float a00 = a.m00;
        final float a01 = a.m01;
        final float a02 = a.m02;
        final float a10 = a.m10;
        final float a11 = a.m11;
        final float a12 = a.m12;
        final float a20 = a.m20;
        final float a21 = a.m21;
        final float a22 = a.m22;
        float b0 = cz*cy;
        float b1 = sz*cy;
        float b2 = -sy;
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

    public AffineTransform3f rotateRPY(float r, float p, float y) {
        return rotateRPY(this, r, p, y);
    }

    public AffineTransform3f setScale(float x, float y, float z) {
        m10 = m20 = m01 = m21 = m02 = m12 = m03 = m13 = m23 = 0.0f;
        m00 = x;
        m11 = y;
        m22 = z;
        return this;
    }

    public AffineTransform3f scale(AffineTransform3f a, float x, float y, float z) {
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

    public AffineTransform3f scale(float x, float y, float z) {
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
    public AffineTransform3f translate(AffineTransform3f a, float x, float y, float z) {
        final float a00 = a.m00;
        final float a10 = a.m10;
        final float a20 = a.m20;
        final float a01 = a.m01;
        final float a11 = a.m11;
        final float a21 = a.m21;
        final float a02 = a.m02;
        final float a12 = a.m12;
        final float a22 = a.m22;

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

    public AffineTransform3f translate(AffineTransform3f a, Vec3f t) {
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
    public AffineTransform3f translate(float x, float y, float z) {
        // inlined translate(this, x, y, z) call:
        this.m03 += m00*x + m01*y + m02*z;
        this.m13 += m10*x + m11*y + m12*z;
        this.m23 += m20*x + m21*y + m22*z;
        return this;
    }

    public AffineTransform3f translate(Vec3f t) {
        return translate(t.x, t.y, t.z);
    }


    /**
     * Sets the rotation to be a random rotation.
     * Leaves the translation untouched.
     *
     * @param rand source of randomness
     * @return this
     */
    public AffineTransform3f randomRotation(Random rand) {
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

//        float x0 = rand.nextFloat();
//        float r1 = (float)Math.sqrt(1.0f - x0);
//        float r2 = (float)Math.sqrt(x0);
//        float t1 = rand.nextFloat()*((float)Math.PI*2);
//        float t2 = rand.nextFloat()*((float)Math.PI*2);
//        float c1 = (float)Math.cos(t1);
//        float s1 = (float)Math.sin(t1);
//        float c2 = (float)Math.cos(t2);
//        float s2 = (float)Math.sin(t2);
//        float qw = s1 * r1;
//        float qx = c1 * r1;
//        float qy = s2 * r2;
//        float qz = c2 * r2;
//        float w2 = qw*2;
//        float wx = w2*qx;
//        float wy = w2*qy;
//        float wz = w2*qz;
//        float x2 = qx*2;
//        float xx = x2*qx;
//        float xy = x2*qy;
//        float xz = x2*qz;
//        float y2 = qy*2;
//        float yy = y2*qy;
//        float yz = y2*qz;
//        float zz = 2*qz*qz;
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
        float r2r2 = rand.nextFloat()*2;
        float r1r1 = (2.0f - r2r2);
        float r1r2 = (float)Math.sqrt(r1r1*r2r2);
        float t1 = rand.nextFloat()*((float)Math.PI*2);
        float t2 = rand.nextFloat()*((float)Math.PI*2);
        float c1 = (float)Math.cos(t1);
        float s1 = (float)Math.sin(t1);
        float c2 = (float)Math.cos(t2);
        float s2 = (float)Math.sin(t2);
        float c1s2 = c1*s2;
        float c2s1 = c2*s1;
        float c1c2 = c1*c2;
        float s1s2 = s1*s2;
        float r1r1c1c1_1 = 1 - r1r1*c1*c1;
        float r2r2c2s2 = r2r2*c2*s2;
        float r1r1c2s1 = r1r1*c1*s1;
        float r2r2s2s2 = r2r2*s2*s2;

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
     * @see #inverse(AffineTransform3f)
     * @return this
     */
    public AffineTransform3f inverse() {
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
    public AffineTransform3f inverse(AffineTransform3f m) {
        // See http://stackoverflow.com/questions/1148309/inverting-a-4x4-matrix
        // This code follows that algorithm, though is significantly reduced based
        // upon the bottom row begin fixed at [0,0,0,1]

        float a11 = m.m11;
        float a21 = m.m21;
        float a22 = m.m22;
        float a12 = m.m12;
        float inv00 = a11*a22 - a12*a21;
        float a01 = m.m01;
        float a02 = m.m02;
        float inv01 = a02*a21 - a01*a22;
        float inv02 = a01*a12 - a02*a11;
        float a00 = m.m00;
        float a10 = m.m10;
        float a20 = m.m20;
        float det = a00*inv00 + a10*inv01 + a20*inv02;

        if (det == 0.0f)
            throw new ArithmeticException("transform not invertible");

        det = 1.0f/det;

        this.m00 = inv00*det;
        this.m01 = inv01*det;
        this.m02 = inv02*det;
        this.m10 = (a12*a20 - a10*a22)*det;
        this.m11 = (a00*a22 - a02*a20)*det;
        this.m12 = (a02*a10 - a00*a12)*det;
        this.m20 = (a10*a21 - a11*a20)*det;
        this.m21 = (a01*a20 - a00*a21)*det;
        this.m22 = (a00*a11 - a01*a10)*det;

        float tx = -m.m03;
        float ty = m.m13;
        float tz = m.m23;
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
     * 1.0f, however the determinant is not a sufficient check for special orthogonality.
     *
     * @return the determinant of the transform matrix
     */
    public float determinant() {
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

        return Math.abs(1 - r00) < 1e-9f &&
               Math.abs(1 - r11) < 1e-9f &&
               Math.abs(1 - r22) < 1e-9f &&
               Math.abs(r01) < 1e-9f &&
               Math.abs(r02) < 1e-9f &&
               Math.abs(r12) < 1e-9f;
    }

    /**
     * Checks if the rotation component of this transform is special orthogonal.
     * A special orthogonal matrix is a proper rotation (no sheers, scales, etc...),
     * and the rotation can be inverted by a transpose.  A matrix is special
     * orthogonal if its is orthogonal and its determinant is 1 (the only other
     * possibility for an orthogonal matrix is -1).
     *
     * <p>In general, the rotation component of an AffineTransform3f instance will
     * remain special orthogonal as long as only rotations (and translations) are
     * applied to the transform.  However numeric precision may cause the matrix
     * to drift out of orthogonality and require fixing.</p>
     *
     * <p>When an AffineTransform3f instance has a special orthogonal rotation
     * component, one can use the fastRigidInverse methods.</p>
     *
     * @return true if the rotation is special orthogonal
     */
    public boolean isRotationSpecialOrthogonal() {
        return Math.abs(1.0f - determinant()) < 1e-9f && isRotationOrthogonal();
    }


    /**
     * Computes an inverse transform with the assumption that
     * the argument has only rotations and translations.
     *
     * @param m the transform to invert
     * @return {@code this}
     */
    public AffineTransform3f fastRigidInverse(AffineTransform3f m) {
        float tx = -m.m03;
        float ty = m.m13;
        float tz = m.m23;
        // in case this == m, store these in temporaries
        float tmp10 = m.m01;
        float tmp20 = m.m02;
        float tmp21 = m.m12;
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
     * @see #fastRigidInverse(AffineTransform3f)
     * @return {@code this}
     */
    public AffineTransform3f fastRigidInverse() {
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
    public AffineTransform3f randomTranslation(
        Random rand,
        float xMin, float xMax,
        float yMin, float yMax,
        float zMin, float zMax)
    {
        this.m03 = rand.nextFloat() * (xMax - xMin) + xMin;
        this.m13 = rand.nextFloat() * (yMax - yMin) + yMin;
        this.m23 = rand.nextFloat() * (zMax - zMin) + zMin;
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
    public AffineTransform3f random(
        Random rand,
        float xMin, float xMax,
        float yMin, float yMax,
        float zMin, float zMax)
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
    public AffineTransform3f rotationTransposeTimes(AffineTransform3f t1, AffineTransform3f t2) {
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
    public AffineTransform3f generateCoordinateSystem(Vec3f w) {
        m00 = w.x;
        m10 = w.y;
        m20 = w.z;

        if (Math.abs(w.x) >= Math.abs(w.y)) {
            float invLen = 1.0f/(float)Math.sqrt(m00*m00 + m20*m20);
            m01 = -m20*invLen;
            m11 = 0;
            m21 = m00*invLen;
            m02 = m10*m21;
            m12 = m20*m01 - m00*m21;
            m22 = -m10*m01;
        } else {
            float inv_length = 1.0f/(float)Math.sqrt(m10*m10 + m20*m20);
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

    public AffineTransform3f setAxes(Vec3f u, Vec3f v, Vec3f w) {
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

    public AffineTransform3f makeRotationOrthogonal() {
        return makeRotationOrthogonal(this);
    }

    public AffineTransform3f makeRotationOrthogonal(AffineTransform3f t) {
        float ux = t.m00;
        float uy = t.m10;
        float uz = t.m20;
        float vx = t.m01;
        float vy = t.m11;
        float vz = t.m21;

        // compute w = u x v
        float wx = uy * vz - uz * vy;
        float wy = uz * vx - ux * vz;
        float wz = ux * vy - uy * vx;

        // normalize w
        float dw = (float)Math.sqrt(wx*wx + wy*wy + wz*wz);
        wx /= dw;
        wy /= dw;
        wz /= dw;

        // normalize v
        float dv = (float)Math.sqrt(vx*vx + vy*vy + vz*vz);
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

    public float dotColumnX(Vec3f v) {
        return m00 * v.x + m10 * v.y + m20 * v.z;
    }

    public float dotColumnY(Vec3f v) {
        return m01 * v.x + m11 * v.y + m21 * v.z;
    }

    public float dotColumnZ(Vec3f v) {
        return m02 * v.x + m12 * v.y + m22 * v.z;
    }

    public float dotColumn(int axis, Vec3f v) {
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
     * {@code +0.0ff == -0.0ff} and {@code NaN != NaN}, then use this
     * method with {@code epsilon = 0.0ff}.  (Contrast with {@link #equals(Object)},
     * which treats {@code +0.0ff != -0.0ff} and {@code NaN == NaN}).
     *
     * @param that the other transform to compare
     * @param epsilon the error range (must be >= 0)
     * @return true if each element is within epsilon of the other.
     */
    public boolean equals(AffineTransform3f that, float epsilon) {
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
    public AffineTransform3f clone() {
        try {
            return (AffineTransform3f)super.clone();
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
        int h =   Float.floatToIntBits(m00);
        h = h*31 + Float.floatToIntBits(m01);
        h = h*31 + Float.floatToIntBits(m02);
        h = h*31 + Float.floatToIntBits(m03);
        h = h*31 + Float.floatToIntBits(m10);
        h = h*31 + Float.floatToIntBits(m11);
        h = h*31 + Float.floatToIntBits(m12);
        h = h*31 + Float.floatToIntBits(m13);
        h = h*31 + Float.floatToIntBits(m20);
        h = h*31 + Float.floatToIntBits(m21);
        h = h*31 + Float.floatToIntBits(m22);
        h = h*31 + Float.floatToIntBits(m23);
        return h ^ HASH_BASE;
    }

    /**
     * Equals method compatible with hashing.  Avoid using this method for comparing
     * two transforms, as it is likely not what you want.  In order to be compatible
     * with hashing contract, the coefficients are compared by their bit values, this
     * means that unlike '==', that NaN's are considered equivalent, and 0.0f is not
     * the same as -0.0f.
     *
     * @param o the object to compare
     * @return true if objects are bit-wise equals.
     */
    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (o == null || o.getClass() != getClass()) return false;
        final AffineTransform3f that = (AffineTransform3f)o;
        return Float.floatToIntBits(this.m00) == Float.floatToIntBits(that.m00)
            && Float.floatToIntBits(this.m01) == Float.floatToIntBits(that.m01)
            && Float.floatToIntBits(this.m02) == Float.floatToIntBits(that.m02)
            && Float.floatToIntBits(this.m03) == Float.floatToIntBits(that.m03)
            && Float.floatToIntBits(this.m10) == Float.floatToIntBits(that.m10)
            && Float.floatToIntBits(this.m11) == Float.floatToIntBits(that.m11)
            && Float.floatToIntBits(this.m12) == Float.floatToIntBits(that.m12)
            && Float.floatToIntBits(this.m13) == Float.floatToIntBits(that.m13)
            && Float.floatToIntBits(this.m20) == Float.floatToIntBits(that.m20)
            && Float.floatToIntBits(this.m21) == Float.floatToIntBits(that.m21)
            && Float.floatToIntBits(this.m22) == Float.floatToIntBits(that.m22)
            && Float.floatToIntBits(this.m23) == Float.floatToIntBits(that.m23);
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
