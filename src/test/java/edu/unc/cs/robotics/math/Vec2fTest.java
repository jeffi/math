// Automatically generated from Vec2dTest.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/4/16.
 */
public class Vec2fTest extends TestCase {
    public void testCtor() throws Exception {
        Vec2f v = new Vec2f();
        assertEquals(0.0f, v.x, 0.0f);
        assertEquals(0.0f, v.y, 0.0f);

        Vec2f u = new Vec2f(1.2f, 3.4f);
        assertEquals(1.2f, u.x, 0.0f);
        assertEquals(3.4f, u.y, 0.0f);
    }
}