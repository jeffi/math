package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/8/16.
 */
public class MathExTest extends TestCase {
    public void testCeilMultiple() throws Exception {
        assertEquals(0.0, MathEx.ceilMultiple(0.0, 10.0), 0.0);
        assertEquals(10.0, MathEx.ceilMultiple(0.000000000001, 10.0), 0.0);
    }
}