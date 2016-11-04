package edu.unc.cs.robotics.math;

/**
 * Created by jeffi on 11/1/16.
 */
public class MatrixIndexOutOfBoundsException extends IndexOutOfBoundsException {
    public MatrixIndexOutOfBoundsException(Matrix m, int r, int c) {
        super(String.format("index (%d, %d) into %dx%d matrix",
            r, c, m.rows(), m.columns()));
    }
}
