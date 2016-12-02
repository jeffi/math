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

    private double _epsilon = 1e-9;
    private int _iterationLimit = 50;

    public NewtonSolver epsilon(double epsilon) {
        _epsilon = epsilon;
        return this;
    }

    public NewtonSolver iterationLimit(int limit) {
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
            if (Math.abs((x1 = x0 - y[0]/y[1]) - x0) < _epsilon) {
                return x1;
            }
        }
        throw new ConvergenceException("iteration limit reached", x1);
    }
}
