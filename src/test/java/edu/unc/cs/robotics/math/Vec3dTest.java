package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/4/16.
 */
public class Vec3dTest extends TestCase {

    private static void assertVec3d(Vec3d n, double x, double y, double z, double delta) {
        assertEquals(x, n.x, delta);
        assertEquals(y, n.y, delta);
        assertEquals(z, n.z, delta);
    }

    /**
     * Test normalizing the unit axis vectors
     */
    public void testNormalizeSimple() throws Exception {
        Vec3d n = new Vec3d();
        assertTrue(n.setCoeffs(0,0,1).normalize());
        assertVec3d(n, 0.0, 0.0, 1f, 0.0);
        assertTrue(n.setCoeffs(0,1,0).normalize());
        assertVec3d(n, 0.0, 1f, 0.0, 0.0);
        assertTrue(n.setCoeffs(1,0,0).normalize());
        assertVec3d(n, 1f, 0.0, 0.0, 0.0);
    }

    /**
     * Try normalizing (1,1,1) at different scales
     */
    public void testNormalizeScale() throws Exception {
        Vec3d n = new Vec3d();
        double r3 = Math.sqrt(1.0/3.0);
        assertTrue(n.setCoeffs(1,1,1).normalize());
        assertVec3d(n, r3, r3, r3, 1e-10);
        assertTrue(n.setCoeffs(-1,-1,-1).normalize());
        assertVec3d(n, -r3, -r3, -r3, 1e-10);
        assertTrue(n.setCoeffs(1e3,-1e3,1e3).normalize());
        assertVec3d(n, r3, -r3, r3, 1e-10);
        assertTrue(n.setCoeffs(1e-3,1e-3,-1e-3).normalize());
        assertVec3d(n, r3, r3, -r3, 1e-10);
    }

    /**
     * Make sure that we are getting the axes correct
     */
    public void testNormalizeCorrectAxes() throws Exception {
        Vec3d n = new Vec3d(3,5,7);
        double d = Math.sqrt(3*3 + 5*5 + 7*7);
        assertTrue(n.normalize());
        assertVec3d(n, 3/d, 5/d, 7/d, 0.0);
    }

    /**
     * Normalizing a zero vector returns false, and
     * does not change the vector.
     */
    public void testNormalizeZero() throws Exception {
        Vec3d n = new Vec3d();
        assertFalse(n.normalize());
        assertVec3d(n, 0.0, 0.0, 0.0, 0.0);
    }

    /**
     * Normalizing a vector with a NaN returns false,
     * and does not change the vector.
     */
    public void testNormalizeWithNaN() throws Exception {
        Vec3d n = new Vec3d(0, 1, Double.NaN);
        assertFalse(n.normalize());
        assertEquals(0.0, n.x, 0.0);
        assertEquals(1f, n.y, 0.0);
        assertEquals(Double.doubleToLongBits(Double.NaN), Double.doubleToLongBits(n.z));
    }

    /**
     * When normalizing a really small vector, we encounter the
     * problem that the squaredNorm() underflows and results in
     * a 0.  Special case handling in normalize() fixes this
     * issue (though there's nothing we can really do about
     * squaredNorm()).
     */
    public void testNormalizeReallySmall() throws Exception {
        Vec3d n = new Vec3d(1e-30, 0, 0);
        assertTrue(n.normalize());
        assertVec3d(n, 1f, 0.0, 0.0, 0.0);
        double d = Math.sqrt(3*3 + 5*5 + 7*7);
        assertTrue(n.setCoeffs(3e-30, 5e-30, 7e-30).normalize());
        assertVec3d(n, 3/d, 5/d, 7/d, 1e-10);
    }

    public void testNormalizeReallyLarge() throws Exception {
        Vec3d n = new Vec3d(3e+30, 5e+30, 7e+30);
        assertTrue(n.normalize());
        double d = Math.sqrt(3*3 + 5*5 + 7*7);
        assertVec3d(n, 3/d, 5/d, 7/d, 1e-7f);
    }

    public void testTransform() throws Exception {
        AffineTransform3d t = new AffineTransform3d();
        t.translate(2, 3, 5);
        t.rotateX(Math.PI / 2);
        Vec3d v = new Vec3d(7, 11, 13);
        Vec3d r = new Vec3d();

        r.transform(t, v);
        assertVec3d(r, 2+7, 3-13, 5+11, 1e-9);
        assertVec3d(v, 7, 11, 13, 0);
    }

    public void testInverseTransform() throws Exception {
        AffineTransform3d t = new AffineTransform3d();
        t.translate(2, 3, 5);
        t.rotateX(Math.PI / 2);
        Vec3d v = new Vec3d(2+7, 3-13, 5+11);
        Vec3d r = new Vec3d();
        r.inverseTransform(t, v);
        assertVec3d(r, 7, 11, 13, 1e-9);
        assertVec3d(v, 2+7, 3-13, 5+11, 1e-9);
    }

    public void testFastRigidInverseTransform() throws Exception {
        AffineTransform3d t = new AffineTransform3d();
        t.translate(2, 3, 5);
        t.rotateX(Math.PI / 2);
        Vec3d v = new Vec3d(2+7, 3-13, 5+11);
        Vec3d r = new Vec3d();
        r.fastRigidInverseTransform(t, v);
        assertVec3d(r, 7, 11, 13, 1e-9);
        assertVec3d(v, 2+7, 3-13, 5+11, 1e-9);
    }
}