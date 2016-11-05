package edu.unc.cs.robotics.math;

/**
 * Thrown when a matrix index pair (r,c) is out of bounds
 */
public class MatrixIndexOutOfBoundsException extends IndexOutOfBoundsException {

    private static final long serialVersionUID = -2395585001320676235L;

    public MatrixIndexOutOfBoundsException(Matrix m, int r, int c) {
        super(String.format("index (%d, %d) into %dx%d matrix",
            r, c, m.rows(), m.columns()));
    }
}
