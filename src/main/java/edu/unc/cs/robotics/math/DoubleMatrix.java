package edu.unc.cs.robotics.math;

/**
 * Typed extension of Matrix.
 */
public interface DoubleMatrix extends Matrix {

    /**
     * Gets the coefficient at the specified row and column.
     * Coefficients are indexed 0-based.
     *
     * @param r the row
     * @param c the column
     * @return the coefficient at (r,c)
     * @throws MatrixIndexOutOfBoundsException if (r,c) is out of bounds.
     */
    double getCoeff(int r, int c);

    /**
     * Sets the value of the coefficient at the specified row and column.
     * Coefficients are indexed 0-based.
     *
     * @param r the row
     * @param c the column
     * @param value the value to set at (r,c)
     * @throws MatrixIndexOutOfBoundsException if (r,c) is out of bounds.
     */
    void setCoeff(int r, int c, double value);
}
