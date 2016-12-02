package edu.unc.cs.robotics.math;

import java.util.concurrent.atomic.AtomicInteger;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/22/16.
 */
public class NewtonSolverTest extends TestCase {
    public void testSqrt() throws Exception {
        // compute sqrt(612)
        // thus, find root of x*x = 612
        // derivative is 2*x
        NewtonSolver solver = new NewtonSolver();
        double root = solver.solve(10, (y, x) -> {
            y[0] = x*x - 612;
            y[1] = 2*x;
        });
        assertEquals(Math.sqrt(612), root, 1e-9);
    }

    public void testCosXeqX3() throws Exception {
        // compute cos(x) = x^3
        NewtonSolver solver = new NewtonSolver();
        double root = solver.solve(0.5, (y, x) -> {
            y[0] = Math.cos(x) - x*x*x;
            y[1] = -Math.sin(x) - 3*x*x;
        });
        assertEquals(0.0, Math.cos(root) - Math.pow(root,3), 1e-9);
    }

    public void testNonConvergence() throws Exception {
        AtomicInteger evaluations = new AtomicInteger(0);
        try {
            // This function has a root, however the Newton method
            // oscillates between 0 and 1.  This test checks that
            // the solver will terminate with an exception rather
            // than infinite loop.
            NewtonSolver solver = new NewtonSolver()
                .iterationLimit(25);
            double root = solver.solve(0.0, (y, x) -> {
                y[0] = x*x*x - 2*x + 2;
                y[1] = 3*x*x - 2;
                evaluations.incrementAndGet();
            });
            fail("converged to "+root+", wasn't expecting that");
        } catch (ConvergenceException ex) {
            // expected
            assertEquals(25, evaluations.get());
        }
    }

    public void testBadEval() throws Exception {
        try {
            NewtonSolver solver = new NewtonSolver();
            solver.solve(0.0, (y, x) -> {
                y[0] = 1;
                y[1] = Double.NaN;
            });
            fail("expected an exception due to the NaN");
        } catch (IllegalStateException ex) {
            // expected
        }

        try {
            NewtonSolver solver = new NewtonSolver();
            solver.solve(0.0, (y, x) -> {
                y[0] = Double.NaN;
                y[1] = 1;
            });
            fail("expected an exception due to the NaN");
        } catch (IllegalStateException ex) {
            // expected
        }

    }
}