package edu.unc.cs.robotics.math;

/**
 * A univariate Newton's method root finder.
 */
public class NewtonSolver {
    @FunctionalInterface
    public interface Fn {
        /**
         * Compute the function and its derivatives and populate the values in x, s.t.,
         * y[0] is the value, and y[1] is the first derivative.
         *
         * @param y the holder for the values.
         * @param x the argument to the function.
         */
        void compute(double[] y, double x);
    }

    private double _tolerance = 1e-9;
    private int _iterationLimit = 50;

    /**
     * Specifies the difference tolerance.  Specify this as 1e-# where # is the
     * number of digits of accuracy desired.  The default is 1e-9 for 9 digits
     * of accuracy.
     *
     * @param epsilon the tolerance required in the result
     * @return {@code this}
     * @throws IllegalArgumentException if argument is less than 0 or greater than 0.1
     */
    public NewtonSolver tolerance(double epsilon) {
        if (epsilon < 0.0 || epsilon >= 1e-1) {
            throw new IllegalArgumentException("tolerance must be between 0 and 0.1");
        }
        _tolerance = epsilon;
        return this;
    }

    /**
     * Specifies the tolerance.
     *
     * @param epsilon the tolerance
     * @return {@code this}
     * @deprecated use {@link #tolerance(double)} instead.
     */
    @Deprecated
    public NewtonSolver epsilon(double epsilon) {
        return tolerance(epsilon);
    }

    /**
     * Specifies the number of iterations to compute before giving up and throwing
     * a {@link ConvergenceException}.
     *
     * @param limit the iteration limit.
     * @return {@code this}
     * @throws IllegalArgumentException if limit is not at least 1
     */
    public NewtonSolver iterationLimit(int limit) {
        if (limit < 1) {
            throw new IllegalArgumentException("iteration limit must be at least 1");
        }
        _iterationLimit = limit;
        return this;
    }

    public double solve(double initialGuess, Fn fn) {
        double x1 = initialGuess, x0;
        double[] y = new double[2];
        for (int i=0, limit=_iterationLimit ; i<limit ; ++i) {
            fn.compute(y, x0 = x1);
            if (Double.isNaN(y[0])) {
                throw new IllegalStateException("value at "+x0+" is NaN");
            }
            if (Double.isNaN(y[1])) {
                throw new IllegalStateException("derivative at "+x0+" is NaN");
            }
            if (y[1] == 0.0) {
                throw new ConvergenceException("derivative at "+x0+" is 0.0", x0);
            }
            x1 = x0 - y[0]/y[1];
            if (Math.abs(x1 - x0) < Math.abs(x1)*_tolerance) {
                return x1;
            }
        }
        throw new ConvergenceException("iteration limit reached", x1);
    }
}
