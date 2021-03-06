package edu.unc.cs.robotics.math;

import java.util.Random;

import junit.framework.TestCase;

/**
 * Created by jeffi on 3/9/16.
 */
public class AffineTransform3dTest extends TestCase {

    static final double EPSILON = 1e-9;
    static final double DEG30 = (double)(30*Math.PI/180);

    static void assertVec3d(double x, double y, double z, Vec3d actual) {
        assertEquals("x", x, actual.x, 1e-6);
        assertEquals("y", y, actual.y, 1e-6);
        assertEquals("z", z, actual.z, 1e-6);
    }

    static void assertAffineTransform3d(
        double m00, double m01, double m02, double m03,
        double m10, double m11, double m12, double m13,
        double m20, double m21, double m22, double m23,
        AffineTransform3d actual)
    {
        final double delta = 1e-6;

        assertEquals("m00", m00, actual.m00, delta);
        assertEquals("m01", m01, actual.m01, delta);
        assertEquals("m02", m02, actual.m02, delta);
        assertEquals("m03", m03, actual.m03, delta);

        assertEquals("m10", m10, actual.m10, delta);
        assertEquals("m11", m11, actual.m11, delta);
        assertEquals("m12", m12, actual.m12, delta);
        assertEquals("m13", m13, actual.m13, delta);

        assertEquals("m20", m20, actual.m20, delta);
        assertEquals("m21", m21, actual.m21, delta);
        assertEquals("m22", m22, actual.m22, delta);
        assertEquals("m23", m23, actual.m23, delta);
    }

    private static void assertValidRotation(AffineTransform3d actual) {
        Matrix3d m = new Matrix3d().fromTransform(actual);
//        System.out.println(m.determinant());
//        System.out.println(m.isOrthogonal());
        assertTrue(m.isSpecialOrthogonal());
    }

    private static void assertAffineTransform3d(AffineTransform3d expected, AffineTransform3d actual) {
        assertAffineTransform3d(
            expected.m00, expected.m01, expected.m02, expected.m03,
            expected.m10, expected.m11, expected.m12, expected.m13,
            expected.m20, expected.m21, expected.m22, expected.m23,
            actual);
    }

    static void assertEquals(
        AffineTransform3d expected,
        AffineTransform3d actual)
    {
        assertAffineTransform3d(
            expected.m00, expected.m01, expected.m02, expected.m03,
            expected.m10, expected.m11, expected.m12, expected.m13,
            expected.m20, expected.m21, expected.m22, expected.m23,
            actual);
    }

    /**
     * Creates and returns a transform where each value is labelled
     * according to its position.  The result is not a valid rigid
     * transform.
     *
     * @return a transform
     */
    static AffineTransform3d createLabelledTransform() {
        return new AffineTransform3d(
            2.00, 0.01, 0.02, 0.03,
            0.10, 0.11, 0.12, 0.13,
            0.20, 0.21, 0.22, 0.23);
    }

    public void testCreateIdentity() throws Exception {
        final AffineTransform3d t = new AffineTransform3d();
        assertAffineTransform3d(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            t);
    }

    public void testCreate() throws Exception {
        assertAffineTransform3d(
            2.00, 0.01, 0.02, 0.03,
            0.10, 0.11, 0.12, 0.13,
            0.20, 0.21, 0.22, 0.23,
            createLabelledTransform());
    }

    public void testSetIdentity() throws Exception {
        final AffineTransform3d t = createLabelledTransform();
        assertAffineTransform3d(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            t.identity());
    }

    public void testSetFromOther() throws Exception {
        final AffineTransform3d t = new AffineTransform3d();
        final AffineTransform3d o = createLabelledTransform();
        t.set(o);
        assertAffineTransform3d(
            2.00, 0.01, 0.02, 0.03,
            0.10, 0.11, 0.12, 0.13,
            0.20, 0.21, 0.22, 0.23,
            t);
    }

    public void testEquals() throws Exception {
        assertEquals(createLabelledTransform(), createLabelledTransform());
        final AffineTransform3d a = new AffineTransform3d();
        final AffineTransform3d b = new AffineTransform3d();
        assertTrue(a.equals(b));
        a.m00 = Double.NaN;
        assertFalse(a.equals(b));
        b.m00 = Double.NaN;
        assertTrue(a.equals(b));
        a.m01 = -0.0;
        b.m01 = +0.0;
        assertFalse(a.equals(b));
        b.m01 = -0.0;
        assertTrue(a.equals(b));
    }

    public void testGet() throws Exception {
        final AffineTransform3d t = createLabelledTransform();
        assertEquals(2.00, t.getCoeff(0, 0), 0.0);
        assertEquals(0.01, t.getCoeff(0, 1), 0.0);
        assertEquals(0.02, t.getCoeff(0, 2), 0.0);
        assertEquals(0.03, t.getCoeff(0, 3), 0.0);

        assertEquals(0.10, t.getCoeff(1, 0), 0.0);
        assertEquals(0.11, t.getCoeff(1, 1), 0.0);
        assertEquals(0.12, t.getCoeff(1, 2), 0.0);
        assertEquals(0.13, t.getCoeff(1, 3), 0.0);

        assertEquals(0.20, t.getCoeff(2, 0), 0.0);
        assertEquals(0.21, t.getCoeff(2, 1), 0.0);
        assertEquals(0.22, t.getCoeff(2, 2), 0.0);
        assertEquals(0.23, t.getCoeff(2, 3), 0.0);
    }

    public void testSet() throws Exception {
        final AffineTransform3d t = new AffineTransform3d();

        t.setCoeff(0, 0, 100);
        t.setCoeff(0, 1, 101);
        t.setCoeff(0, 2, 102);
        t.setCoeff(0, 3, 103);

        t.setCoeff(1, 0, 110);
        t.setCoeff(1, 1, 111);
        t.setCoeff(1, 2, 112);
        t.setCoeff(1, 3, 113);

        t.setCoeff(2, 0, 120);
        t.setCoeff(2, 1, 121);
        t.setCoeff(2, 2, 122);
        t.setCoeff(2, 3, 123);

        assertAffineTransform3d(
            100, 101, 102, 103,
            110, 111, 112, 113,
            120, 121, 122, 123,
            t);
    }

    public void testHashCode() throws Exception {
        // create a zero matrix by setting the diagonal to 0
        final AffineTransform3d t = new AffineTransform3d();
        final AffineTransform3d o = new AffineTransform3d();
        t.m00 = t.m11 = t.m22 = 0.0;

        o.set(t);

        // hash base is used to make sure that we don't get
        // the same hash code for things like 0 objects of
        // different types
        final int hashBase = AffineTransform3d.class.getName().hashCode();

        // If all values are 0, then the hashcode will be 0.
        assertEquals(0 ^ hashBase, t.hashCode());
        assertTrue(o.equals(t));

        // -0.0 is 0.0 with the sign bit set, this the bit
        // pattern will be 1000....0.  By setting the last
        // entry (m23) we should get that value out.  This
        // tests that +0.0 and -0.0 are treated differently.
        t.m23 = -0.0;
        assertEquals(0x80000000 ^ hashBase, t.hashCode());
        assertFalse(o.equals(t));

        // On the other hand, NaN != NaN, but for hashing
        // purposes hash(NaN) == hash(NaN).
        t.m23 = o.m23 = Double.NaN;
        int prev = new Double(Double.NaN).hashCode();
//        int prev = 0x7ff80000; // Double.doubleToLongBits(double.NaN) shifted and xor'd;
        assertEquals(prev ^ hashBase, t.hashCode());
        assertTrue(o.equals(t));

        // Now we loop through each entry and check that
        // changing a single value changes the hash code.
        for (int r = 0; r < 3; ++r) {
            for (int c = 0; c < 4; ++c) {
                o.set(t);
                assertTrue(o.equals(t));
                t.setCoeff(r, c, 1);
                final int curr = t.hashCode();
                assertFalse(prev == curr);
                assertFalse(o.equals(t));
                prev = curr;
            }
        }
    }

    public void testRotateX() throws Exception {
        // rotate 30 degrees around x-axis
        for (int i=0 ; i<4 ; ++i) {
            final AffineTransform3d t = new AffineTransform3d();
            switch (i) {
            case 0:
                t.rotationX(DEG30);
                break;
            case 1:
                t.rotationX(DEG30/2);
                t.rotateX(DEG30/2);
                break;
            case 2:
                t.rotation(1, 0, 0, DEG30);
                break;
            case 3:
                t.rotationRPY(DEG30, 0, 0);
                break;
            }

            assertVec3d(1.0, 0.8660254, 0.5, new Vec3d().transform(t, new Vec3d(1, 1, 0)));
            // y = cos(30) - sin(30) = 0.8660 - 0.5000 = 0.3660
            // z = sin(30) + cos(30) = 0.5000 + 0.8660 = 1.3660
            assertVec3d(1.0, 0.3660254, 1.3660254, new Vec3d().transform(t, new Vec3d(1, 1, 1)));
            assertValidRotation(t);
        }
    }

    public void testRotateY() throws Exception {
        // rotate 30 degrees around the y-axis
        for (int i=0 ; i<4; ++i) {
            final AffineTransform3d t = new AffineTransform3d();
            switch (i) {
            case 0:
                t.rotationY(DEG30);
                break;
            case 1:
                t.rotationY(DEG30/2);
                t.rotateY(DEG30/2);
                break;
            case 2:
                t.rotation(0, 1, 0, DEG30);
                break;
            case 3:
                t.rotationRPY(0, DEG30, 0);
                break;
            }

            assertVec3d(0.8660254, 1.0, -0.5, new Vec3d().transform(t, new Vec3d(1, 1, 0)));
            assertVec3d(1.3660254, 1.0, 0.3660254, new Vec3d().transform(t, new Vec3d(1, 1, 1)));
            assertValidRotation(t);
        }
    }

    public void testRotateZ() throws Exception {
        // rotate 30 degrees around the z-axis
        for (int i=0 ; i<5 ; ++i) {
            final AffineTransform3d t = new AffineTransform3d();
            switch (i) {
            case 0:
                t.rotationZ(DEG30);
                break;
            case 1:
                t.rotationZ(DEG30/2);
                t.rotateZ(DEG30/2);
                break;
            case 2:
                t.rotateZ(new AffineTransform3d().rotationZ(DEG30/2), DEG30/2);
                break;
            case 3:
                t.rotation(0, 0, 1, DEG30);
                break;
            case 4:
                t.rotationRPY(0, 0, DEG30);
                break;
            }

            assertVec3d(0.8660254, 0.5, 1.0, new Vec3d().transform(t, new Vec3d(1, 0, 1)));
            assertVec3d(0.3660254, 1.3660254, 1.0, new Vec3d().transform(t, new Vec3d(1, 1, 1)));
            assertValidRotation(t);
        }
    }

    public void testSetRotate() throws Exception {
        Random rand = new Random(1);
        AffineTransform3d tf = new AffineTransform3d();
        for (int i=0 ; i<1000 ; ++i) {
            tf.rotation(rand.nextDouble()*2-1, rand.nextDouble()*2-1, rand.nextDouble()*2-1, 2*Math.PI*rand.nextDouble());
            assertValidRotation(tf);
        }
    }

    public void testSetRPY() throws Exception {
        AffineTransform3d expected = new AffineTransform3d();
        expected.rotateZ(0.3);
        expected.rotateY(0.2);
        expected.rotateX(0.1);

        AffineTransform3d actual = new AffineTransform3d();
        actual.rotationRPY(0.1, 0.2, 0.3);

        assertEquals(expected, actual);
        assertValidRotation(actual);
    }

    public void testFromQuaternionIdentity() throws Exception {
        Quaternion4d q = new Quaternion4d();
        AffineTransform3d e = new AffineTransform3d();
        AffineTransform3d a = new AffineTransform3d();
        a.rotation(q);
        assertAffineTransform3d(e, a);
    }

    public void testFromQuaternionRPY() throws Exception {
        Quaternion4d q = new Quaternion4d();
        q.fromRPY(1, 2, 3);
        AffineTransform3d e = new AffineTransform3d();
        e.rotationRPY(1, 2, 3);
        AffineTransform3d a = new AffineTransform3d();
        a.rotation(q);
        assertValidRotation(e);
        assertValidRotation(a);
        assertAffineTransform3d(e, a);
    }

    public void testRandomRotation() throws Exception {
        for (int seed=1 ; seed<=100 ; ++seed) {
            Quaternion4d q = new Quaternion4d();
            q.random(new Random(seed));
            assertEquals(1, q.dot(q), EPSILON);
            AffineTransform3d e = new AffineTransform3d();
            e.rotation(q);
            assertValidRotation(e);
            AffineTransform3d a = new AffineTransform3d();
            a.randomRotation(new Random(seed));
            assertValidRotation(a);
            assertEquals(e, a);
        }
    }

    public void testInverse() throws Exception {
        AffineTransform3d T = new AffineTransform3d()
            .rotate(0.1, 2, -3, 1.3)
            .translate(2, -4, 3);
        AffineTransform3d I = new AffineTransform3d().inverse(T);

        assertEquals(T.m00, I.m00, EPSILON);
        assertEquals(I.m01, T.m10, EPSILON);
        assertEquals(I.m02, T.m20, EPSILON);
        assertEquals(I.m10, T.m01, EPSILON);
        assertEquals(I.m11, T.m11, EPSILON);
        assertEquals(I.m12, T.m21, EPSILON);
        assertEquals(I.m20, T.m02, EPSILON);
        assertEquals(I.m21, T.m12, EPSILON);
        assertEquals(I.m22, T.m22, EPSILON);

        assertEquals(-2.0, I.m03, EPSILON);
        assertEquals(4.0, I.m13, EPSILON);
        assertEquals(-3.0, I.m23, EPSILON);
    }

    public void testInverseWithScale() throws Exception {
        AffineTransform3d T = new AffineTransform3d()
            .rotate(0.1, 2, -3, 1.3)
            .translate(2, -4, 3)
            .scale(5,7,11);
        AffineTransform3d I = new AffineTransform3d().inverse(T);

        AffineTransform3d R = T.mul(I);
        AffineTransform3d identity = new AffineTransform3d();
        assertEquals(identity, R);
    }

    public void testFastRigidInverse() throws Exception {
        AffineTransform3d T = new AffineTransform3d()
            .rotate(0.1, 2, -3, 1.3)
            .translate(2, -4, 3);
        AffineTransform3d I = new AffineTransform3d().fastRigidInverse(T);

        assertEquals(T.m00, I.m00, EPSILON);
        assertEquals(I.m01, T.m10, EPSILON);
        assertEquals(I.m02, T.m20, EPSILON);
        assertEquals(I.m10, T.m01, EPSILON);
        assertEquals(I.m11, T.m11, EPSILON);
        assertEquals(I.m12, T.m21, EPSILON);
        assertEquals(I.m20, T.m02, EPSILON);
        assertEquals(I.m21, T.m12, EPSILON);
        assertEquals(I.m22, T.m22, EPSILON);

        assertEquals(-2.0, I.m03, EPSILON);
        assertEquals(4.0, I.m13, EPSILON);
        assertEquals(-3.0, I.m23, EPSILON);
    }

    public void testBenchmark() {
        if (!"true".equals(System.getProperty("benchmark"))) {
            System.out.println("skipping benchmarks, set -Dbenchmark=true to enable");
            return;
        }

        // On MacBook Pro, this benchmarks at ~12.5 ns/mul
        // Bizarrely, the same ops with Eigen::Transform run at approx ~24 ns/mul
        // When the same transform multiply code is reproduced in C++ directly
        // then we see about ~11.9 ns/mul in C/C++ (w/clang++ -O3 -favx -ffma)

        AffineTransform3d a = new AffineTransform3d();
        AffineTransform3d b = new AffineTransform3d();
        AffineTransform3d c = new AffineTransform3d();

        final int N = 1000000;
        double sum = 0.0;
        double count = 0;
        for (int i=0 ; i < 100 ; ++i) {
            c.identity();
            a.rotation(-1, 2, 3, 0.1);
            b.rotation(3, 1, -2, 0.1);

            long startTime = System.nanoTime();
            bench(a, b, c, N);
            long elapsed = System.nanoTime() - startTime;
            double nsPerOp = elapsed/(2.0*N);
            System.out.printf("%.1f ns/mul%n", nsPerOp);
            if (i > 5) {
                sum += nsPerOp;
                count++;
            }
        }
        System.out.println("Average: "+sum/count);
        System.out.println(c);
    }

    private static void bench(AffineTransform3d a, AffineTransform3d b, AffineTransform3d c, int N) {
        for (int i = 0 ; i < N ; ++i) {
            c.mul(c, a);
            c.mul(c, b);
        }
    }
}
