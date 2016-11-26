package edu.unc.cs.robotics.math;

/**
 * Exception thrown with the solver fails to converge.
 */
public class ConvergenceException extends RuntimeException {

    private final double _value;

    public ConvergenceException(double value) {
        _value = value;
    }

    /**
     * Returns the value that did not converge as the solver computed
     * before giving up.  This value may be meaningless.
     *
     * @return the value
     */
    public double getValue() {
        return _value;
    }
}
