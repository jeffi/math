// Automatically generated from AffineTransform3dTest.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import java.util.Random;

import junit.framework.TestCase;

/**
 * Created by jeffi on 3/9/16.
 */
public class AffineTransform3fTest extends TestCase {

    static final float EPSILON = 1e-6f;
    static final float DEG30 = (float)(30*(float)Math.PI/180);

    static void assertVec3f(float x, float y, float z, Vec3f actual) {
        assertEquals("x", x, actual.x, 1e-6f);
        assertEquals("y", y, actual.y, 1e-6f);
        assertEquals("z", z, actual.z, 1e-6f);
    }

    static void assertAffineTransform3f(
        float m00, float m01, float m02, float m03,
        float m10, float m11, float m12, float m13,
        float m20, float m21, float m22, float m23,
        AffineTransform3f actual)
    {
        final float delta = 1e-6f;

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

    private static void assertValidRotation(AffineTransform3f actual) {
        Matrix3f m = new Matrix3f().fromTransform(actual);
//        System.out.println(m.determinant());
//        System.out.println(m.isOrthogonal());
        assertTrue(m.isSpecialOrthogonal());
    }

    private static void assertAffineTransform3f(AffineTransform3f expected, AffineTransform3f actual) {
        assertAffineTransform3f(
            expected.m00, expected.m01, expected.m02, expected.m03,
            expected.m10, expected.m11, expected.m12, expected.m13,
            expected.m20, expected.m21, expected.m22, expected.m23,
            actual);
    }

    static void assertEquals(
        AffineTransform3f expected,
        AffineTransform3f actual)
    {
        assertAffineTransform3f(
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
    static AffineTransform3f createLabelledTransform() {
        return new AffineTransform3f(
            2.00f, 0.01f, 0.02f, 0.03f,
            0.10f, 0.11f, 0.12f, 0.13f,
            0.20f, 0.21f, 0.22f, 0.23f);
    }

    public void testCreateIdentity() throws Exception {
        final AffineTransform3f t = new AffineTransform3f();
        assertAffineTransform3f(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            t);
    }

    public void testCreate() throws Exception {
        assertAffineTransform3f(
            2.00f, 0.01f, 0.02f, 0.03f,
            0.10f, 0.11f, 0.12f, 0.13f,
            0.20f, 0.21f, 0.22f, 0.23f,
            createLabelledTransform());
    }

    public void testSetIdentity() throws Exception {
        final AffineTransform3f t = createLabelledTransform();
        assertAffineTransform3f(
            1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            t.identity());
    }

    public void testSetFromOther() throws Exception {
        final AffineTransform3f t = new AffineTransform3f();
        final AffineTransform3f o = createLabelledTransform();
        t.set(o);
        assertAffineTransform3f(
            2.00f, 0.01f, 0.02f, 0.03f,
            0.10f, 0.11f, 0.12f, 0.13f,
            0.20f, 0.21f, 0.22f, 0.23f,
            t);
    }

    public void testEquals() throws Exception {
        assertEquals(createLabelledTransform(), createLabelledTransform());
        final AffineTransform3f a = new AffineTransform3f();
        final AffineTransform3f b = new AffineTransform3f();
        assertTrue(a.equals(b));
        a.m00 = Float.NaN;
        assertFalse(a.equals(b));
        b.m00 = Float.NaN;
        assertTrue(a.equals(b));
        a.m01 = -0.0f;
        b.m01 = +0.0f;
        assertFalse(a.equals(b));
        b.m01 = -0.0f;
        assertTrue(a.equals(b));
    }

    public void testGet() throws Exception {
        final AffineTransform3f t = createLabelledTransform();
        assertEquals(2.00f, t.getCoeff(0, 0), 0.0f);
        assertEquals(0.01f, t.getCoeff(0, 1), 0.0f);
        assertEquals(0.02f, t.getCoeff(0, 2), 0.0f);
        assertEquals(0.03f, t.getCoeff(0, 3), 0.0f);

        assertEquals(0.10f, t.getCoeff(1, 0), 0.0f);
        assertEquals(0.11f, t.getCoeff(1, 1), 0.0f);
        assertEquals(0.12f, t.getCoeff(1, 2), 0.0f);
        assertEquals(0.13f, t.getCoeff(1, 3), 0.0f);

        assertEquals(0.20f, t.getCoeff(2, 0), 0.0f);
        assertEquals(0.21f, t.getCoeff(2, 1), 0.0f);
        assertEquals(0.22f, t.getCoeff(2, 2), 0.0f);
        assertEquals(0.23f, t.getCoeff(2, 3), 0.0f);
    }

    public void testSet() throws Exception {
        final AffineTransform3f t = new AffineTransform3f();

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

        assertAffineTransform3f(
            100, 101, 102, 103,
            110, 111, 112, 113,
            120, 121, 122, 123,
            t);
    }

    public void testHashCode() throws Exception {
        // create a zero matrix by setting the diagonal to 0
        final AffineTransform3f t = new AffineTransform3f();
        final AffineTransform3f o = new AffineTransform3f();
        t.m00 = t.m11 = t.m22 = 0.0f;

        o.set(t);

        // hash base is used to make sure that we don't get
        // the same hash code for things like 0 objects of
        // different types
        final int hashBase = AffineTransform3f.class.getName().hashCode();

        // If all values are 0, then the hashcode will be 0.
        assertEquals(0 ^ hashBase, t.hashCode());
        assertTrue(o.equals(t));

        // -0.0f is 0.0f with the sign bit set, this the bit
        // pattern will be 1000....0.  By setting the last
        // entry (m23) we should get that value out.  This
        // tests that +0.0f and -0.0f are treated differently.
        t.m23 = -0.0f;
        assertEquals(0x80000000 ^ hashBase, t.hashCode());
        assertFalse(o.equals(t));

        // On the other hand, NaN != NaN, but for hashing
        // purposes hash(NaN) == hash(NaN).
        t.m23 = o.m23 = Float.NaN;
        int prev = new Float(Float.NaN).hashCode();
//        int prev = 0x7ff80000; // Float.floatToIntBits(float.NaN) shifted and xor'd;
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
            final AffineTransform3f t = new AffineTransform3f();
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

            assertVec3f(1.0f, 0.8660254f, 0.5f, new Vec3f().transform(t, new Vec3f(1, 1, 0)));
            // y = cos(30) - sin(30) = 0.8660f - 0.5000f = 0.3660f
            // z = sin(30) + cos(30) = 0.5000f + 0.8660f = 1.3660f
            assertVec3f(1.0f, 0.3660254f, 1.3660254f, new Vec3f().transform(t, new Vec3f(1, 1, 1)));
            assertValidRotation(t);
        }
    }

    public void testRotateY() throws Exception {
        // rotate 30 degrees around the y-axis
        for (int i=0 ; i<4; ++i) {
            final AffineTransform3f t = new AffineTransform3f();
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

            assertVec3f(0.8660254f, 1.0f, -0.5f, new Vec3f().transform(t, new Vec3f(1, 1, 0)));
            assertVec3f(1.3660254f, 1.0f, 0.3660254f, new Vec3f().transform(t, new Vec3f(1, 1, 1)));
            assertValidRotation(t);
        }
    }

    public void testRotateZ() throws Exception {
        // rotate 30 degrees around the z-axis
        for (int i=0 ; i<5 ; ++i) {
            final AffineTransform3f t = new AffineTransform3f();
            switch (i) {
            case 0:
                t.rotationZ(DEG30);
                break;
            case 1:
                t.rotationZ(DEG30/2);
                t.rotateZ(DEG30/2);
                break;
            case 2:
                t.rotateZ(new AffineTransform3f().rotationZ(DEG30/2), DEG30/2);
                break;
            case 3:
                t.rotation(0, 0, 1, DEG30);
                break;
            case 4:
                t.rotationRPY(0, 0, DEG30);
                break;
            }

            assertVec3f(0.8660254f, 0.5f, 1.0f, new Vec3f().transform(t, new Vec3f(1, 0, 1)));
            assertVec3f(0.3660254f, 1.3660254f, 1.0f, new Vec3f().transform(t, new Vec3f(1, 1, 1)));
            assertValidRotation(t);
        }
    }

    public void testSetRotate() throws Exception {
        Random rand = new Random(1);
        AffineTransform3f tf = new AffineTransform3f();
        for (int i=0 ; i<1000 ; ++i) {
            tf.rotation(rand.nextFloat()*2-1, rand.nextFloat()*2-1, rand.nextFloat()*2-1, 2*(float)Math.PI*rand.nextFloat());
            assertValidRotation(tf);
        }
    }

    public void testSetRPY() throws Exception {
        AffineTransform3f expected = new AffineTransform3f();
        expected.rotateZ(0.3f);
        expected.rotateY(0.2f);
        expected.rotateX(0.1f);

        AffineTransform3f actual = new AffineTransform3f();
        actual.rotationRPY(0.1f, 0.2f, 0.3f);

        assertEquals(expected, actual);
        assertValidRotation(actual);
    }

    public void testFromQuaternionIdentity() throws Exception {
        Quaternion4f q = new Quaternion4f();
        AffineTransform3f e = new AffineTransform3f();
        AffineTransform3f a = new AffineTransform3f();
        a.rotation(q);
        assertAffineTransform3f(e, a);
    }

    public void testFromQuaternionRPY() throws Exception {
        Quaternion4f q = new Quaternion4f();
        q.fromRPY(1, 2, 3);
        AffineTransform3f e = new AffineTransform3f();
        e.rotationRPY(1, 2, 3);
        AffineTransform3f a = new AffineTransform3f();
        a.rotation(q);
        assertValidRotation(e);
        assertValidRotation(a);
        assertAffineTransform3f(e, a);
    }

    public void testRandomRotation() throws Exception {
        for (int seed=1 ; seed<=100 ; ++seed) {
            Quaternion4f q = new Quaternion4f();
            q.random(new Random(seed));
            assertEquals(1, q.dot(q), EPSILON);
            AffineTransform3f e = new AffineTransform3f();
            e.rotation(q);
            assertValidRotation(e);
            AffineTransform3f a = new AffineTransform3f();
            a.randomRotation(new Random(seed));
            assertValidRotation(a);
            assertEquals(e, a);
        }
    }

    public void testInverse() throws Exception {
        AffineTransform3f T = new AffineTransform3f()
            .rotate(0.1f, 2, -3, 1.3f)
            .translate(2, -4, 3);
        AffineTransform3f I = new AffineTransform3f().inverse(T);

        assertEquals(T.m00, I.m00, EPSILON);
        assertEquals(I.m01, T.m10, EPSILON);
        assertEquals(I.m02, T.m20, EPSILON);
        assertEquals(I.m10, T.m01, EPSILON);
        assertEquals(I.m11, T.m11, EPSILON);
        assertEquals(I.m12, T.m21, EPSILON);
        assertEquals(I.m20, T.m02, EPSILON);
        assertEquals(I.m21, T.m12, EPSILON);
        assertEquals(I.m22, T.m22, EPSILON);

        assertEquals(-2.0f, I.m03, EPSILON);
        assertEquals(4.0f, I.m13, EPSILON);
        assertEquals(-3.0f, I.m23, EPSILON);
    }

    public void testInverseWithScale() throws Exception {
        AffineTransform3f T = new AffineTransform3f()
            .rotate(0.1f, 2, -3, 1.3f)
            .translate(2, -4, 3)
            .scale(5,7,11);
        AffineTransform3f I = new AffineTransform3f().inverse(T);

        AffineTransform3f R = T.mul(I);
        AffineTransform3f identity = new AffineTransform3f();
        assertEquals(identity, R);
    }

    public void testFastRigidInverse() throws Exception {
        AffineTransform3f T = new AffineTransform3f()
            .rotate(0.1f, 2, -3, 1.3f)
            .translate(2, -4, 3);
        AffineTransform3f I = new AffineTransform3f().fastRigidInverse(T);

        assertEquals(T.m00, I.m00, EPSILON);
        assertEquals(I.m01, T.m10, EPSILON);
        assertEquals(I.m02, T.m20, EPSILON);
        assertEquals(I.m10, T.m01, EPSILON);
        assertEquals(I.m11, T.m11, EPSILON);
        assertEquals(I.m12, T.m21, EPSILON);
        assertEquals(I.m20, T.m02, EPSILON);
        assertEquals(I.m21, T.m12, EPSILON);
        assertEquals(I.m22, T.m22, EPSILON);

        assertEquals(-2.0f, I.m03, EPSILON);
        assertEquals(4.0f, I.m13, EPSILON);
        assertEquals(-3.0f, I.m23, EPSILON);
    }

    public void testBenchmark() {
        if (!"true".equals(System.getProperty("benchmark"))) {
            System.out.println("skipping benchmarks, set -Dbenchmark=true to enable");
            return;
        }

        // On MacBook Pro, this benchmarks at ~12.5f ns/mul
        // Bizarrely, the same ops with Eigen::Transform run at approx ~24 ns/mul
        // When the same transform multiply code is reproduced in C++ directly
        // then we see about ~11.9f ns/mul in C/C++ (w/clang++ -O3 -favx -ffma)

        AffineTransform3f a = new AffineTransform3f();
        AffineTransform3f b = new AffineTransform3f();
        AffineTransform3f c = new AffineTransform3f();

        final int N = 1000000;
        float sum = 0.0f;
        float count = 0;
        for (int i=0 ; i < 100 ; ++i) {
            c.identity();
            a.rotation(-1, 2, 3, 0.1f);
            b.rotation(3, 1, -2, 0.1f);

            long startTime = System.nanoTime();
            bench(a, b, c, N);
            long elapsed = System.nanoTime() - startTime;
            float nsPerOp = elapsed/(2.0f*N);
            System.out.printf("%.1f ns/mul%n", nsPerOp);
            if (i > 5) {
                sum += nsPerOp;
                count++;
            }
        }
        System.out.println("Average: "+sum/count);
        System.out.println(c);
    }

    private static void bench(AffineTransform3f a, AffineTransform3f b, AffineTransform3f c, int N) {
        for (int i = 0 ; i < N ; ++i) {
            c.mul(c, a);
            c.mul(c, b);
        }
    }
}
