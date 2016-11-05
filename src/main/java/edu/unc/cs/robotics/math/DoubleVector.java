package edu.unc.cs.robotics.math;

/**
 * Typed extension of Vector.
 */
public interface DoubleVector extends DoubleMatrix, Vector {
    double getCoeff(int i);
    void setCoeff(int i, double v);

    @Override
    default double getCoeff(int r, int c) {
        if (r != 0) {
            throw new MatrixIndexOutOfBoundsException(this, r, c);
        }
        return getCoeff(c);
    }

    @Override
    default void setCoeff(int r, int c, double value) {
        if (r != 0) {
            throw new MatrixIndexOutOfBoundsException(this, r, c);
        }
        setCoeff(c, value);
    }
}
