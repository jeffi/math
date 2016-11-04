package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/4/16.
 */
public class Vec2dTest extends TestCase {
    public void testCtor() throws Exception {
        Vec2d v = new Vec2d();
        assertEquals(0.0, v.x, 0.0);
        assertEquals(0.0, v.y, 0.0);

        Vec2d u = new Vec2d(1.2, 3.4);
        assertEquals(1.2, u.x, 0.0);
        assertEquals(3.4, u.y, 0.0);
    }
}