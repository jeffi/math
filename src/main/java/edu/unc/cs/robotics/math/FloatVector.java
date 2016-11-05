// Automatically generated from DoubleVector.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

/**
 * Typed extension of Vector.
 */
public interface FloatVector extends FloatMatrix, Vector {
    float getCoeff(int i);
    void setCoeff(int i, float v);

    @Override
    default float getCoeff(int r, int c) {
        if (r != 0) {
            throw new MatrixIndexOutOfBoundsException(this, r, c);
        }
        return getCoeff(c);
    }

    @Override
    default void setCoeff(int r, int c, float value) {
        if (r != 0) {
            throw new MatrixIndexOutOfBoundsException(this, r, c);
        }
        setCoeff(c, value);
    }
}
