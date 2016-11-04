// Automatically generated from Vec3dTest.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/4/16.
 */
public class Vec3fTest extends TestCase {
    
    static final float EPSILON = 1e-6f;

    private static void assertVec3f(Vec3f n, float x, float y, float z, float delta) {
        assertEquals(x, n.x, delta);
        assertEquals(y, n.y, delta);
        assertEquals(z, n.z, delta);
    }

    /**
     * Test normalizing the unit axis vectors
     */
    public void testNormalizeSimple() throws Exception {
        Vec3f n = new Vec3f();
        assertTrue(n.setCoeffs(0,0,1).normalize());
        assertVec3f(n, 0.0f, 0.0f, 1f, 0.0f);
        assertTrue(n.setCoeffs(0,1,0).normalize());
        assertVec3f(n, 0.0f, 1f, 0.0f, 0.0f);
        assertTrue(n.setCoeffs(1,0,0).normalize());
        assertVec3f(n, 1f, 0.0f, 0.0f, 0.0f);
    }

    /**
     * Try normalizing (1,1,1) at different scales
     */
    public void testNormalizeScale() throws Exception {
        Vec3f n = new Vec3f();
        float r3 = (float)Math.sqrt(1.0f/3.0f);
        assertTrue(n.setCoeffs(1,1,1).normalize());
        assertVec3f(n, r3, r3, r3, 1e-10f);
        assertTrue(n.setCoeffs(-1,-1,-1).normalize());
        assertVec3f(n, -r3, -r3, -r3, 1e-10f);
        assertTrue(n.setCoeffs(1e3f,-1e3f,1e3f).normalize());
        assertVec3f(n, r3, -r3, r3, 1e-10f);
        assertTrue(n.setCoeffs(1e-3f,1e-3f,-1e-3f).normalize());
        assertVec3f(n, r3, r3, -r3, 1e-10f);
    }

    /**
     * Make sure that we are getting the axes correct
     */
    public void testNormalizeCorrectAxes() throws Exception {
        Vec3f n = new Vec3f(3,5,7);
        float d = (float)Math.sqrt(3*3 + 5*5 + 7*7);
        assertTrue(n.normalize());
        assertVec3f(n, 3/d, 5/d, 7/d, 0.0f);
    }

    /**
     * Normalizing a zero vector returns false, and
     * does not change the vector.
     */
    public void testNormalizeZero() throws Exception {
        Vec3f n = new Vec3f();
        assertFalse(n.normalize());
        assertVec3f(n, 0.0f, 0.0f, 0.0f, 0.0f);
    }

    /**
     * Normalizing a vector with a NaN returns false,
     * and does not change the vector.
     */
    public void testNormalizeWithNaN() throws Exception {
        Vec3f n = new Vec3f(0, 1, Float.NaN);
        assertFalse(n.normalize());
        assertEquals(0.0f, n.x, 0.0f);
        assertEquals(1f, n.y, 0.0f);
        assertEquals(Float.floatToIntBits(Float.NaN), Float.floatToIntBits(n.z));
    }

    /**
     * When normalizing a really small vector, we encounter the
     * problem that the squaredNorm() underflows and results in
     * a 0.  Special case handling in normalize() fixes this
     * issue (though there's nothing we can really do about
     * squaredNorm()).
     */
    public void testNormalizeReallySmall() throws Exception {
        Vec3f n = new Vec3f(1e-30f, 0, 0);
        assertTrue(n.normalize());
        assertVec3f(n, 1f, 0.0f, 0.0f, 0.0f);
        float d = (float)Math.sqrt(3*3 + 5*5 + 7*7);
        assertTrue(n.setCoeffs(3e-30f, 5e-30f, 7e-30f).normalize());
        assertVec3f(n, 3/d, 5/d, 7/d, 1e-10f);
    }

    public void testNormalizeReallyLarge() throws Exception {
        Vec3f n = new Vec3f(3e+30f, 5e+30f, 7e+30f);
        assertTrue(n.normalize());
        float d = (float)Math.sqrt(3*3 + 5*5 + 7*7);
        assertVec3f(n, 3/d, 5/d, 7/d, 1e-7f);
    }

    public void testTransform() throws Exception {
        AffineTransform3f t = new AffineTransform3f();
        t.translate(2, 3, 5);
        t.rotateX((float)Math.PI / 2);
        Vec3f v = new Vec3f(7, 11, 13);
        Vec3f r = new Vec3f();

        r.transform(t, v);
        assertVec3f(r, 2+7, 3-13, 5+11, EPSILON);
        assertVec3f(v, 7, 11, 13, 0);
    }

    public void testInverseTransform() throws Exception {
        AffineTransform3f t = new AffineTransform3f();
        t.translate(2, 3, 5);
        t.rotateX((float)Math.PI / 2);
        Vec3f v = new Vec3f(2+7, 3-13, 5+11);
        Vec3f r = new Vec3f();
        r.inverseTransform(t, v);
        assertVec3f(r, 7, 11, 13, EPSILON);
        assertVec3f(v, 2+7, 3-13, 5+11, EPSILON);
    }

    public void testFastRigidInverseTransform() throws Exception {
        AffineTransform3f t = new AffineTransform3f();
        t.translate(2, 3, 5);
        t.rotateX((float)Math.PI / 2);
        Vec3f v = new Vec3f(2+7, 3-13, 5+11);
        Vec3f r = new Vec3f();
        r.fastRigidInverseTransform(t, v);
        assertVec3f(r, 7, 11, 13, EPSILON);
        assertVec3f(v, 2+7, 3-13, 5+11, EPSILON);
    }
}