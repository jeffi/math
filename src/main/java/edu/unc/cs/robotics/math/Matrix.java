package edu.unc.cs.robotics.math;

/**
 * Interface of methods common to matrices.
 *
 * @see DoubleMatrix
 * @see FloatMatrix
 * @see Vector
 */
public interface Matrix extends java.io.Serializable {
    /**
     * The number of rows in this matrix
     *
     * @return the number of rows in this matrx
     */
    int rows();

    /**
     * The number of columns in this matrix
     *
     * @return the number of columns in this matrix
     */
    int columns();

    /**
     * The number of coefficients in this matrix (@{code rows()*columns()})
     *
     * @return the number of coefficients in this matrix
     */
    int size();
}
