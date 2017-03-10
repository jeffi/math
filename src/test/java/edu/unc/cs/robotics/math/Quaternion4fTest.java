// Automatically generated from Quaternion4dTest.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import java.util.Random;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/4/16.
 */
public class Quaternion4fTest extends TestCase {
    static final float EPSILON = 1e-6f;

    private static void assertQuaternion(Quaternion4f actual, float ... expected) {
        boolean posEq = true;
        boolean negEq = true;
        float sum = 0;
        for (int i=0 ; i<4 ; ++i) {
            float v = actual.getCoeff(i);
            posEq &= Math.abs(expected[i] - v) < EPSILON;
            negEq &= Math.abs(expected[i] + v) < EPSILON;
            sum += v*v;
        }
        if (!posEq && !negEq) {
            failNotEquals("", new Quaternion4f(expected), actual);
        }
        assertEquals(1.0f, sum, EPSILON);
    }

    public void testCreateIdentity() throws Exception {
        final Quaternion4f q = new Quaternion4f();
        assertEquals(1f, q.w, 0);
        assertEquals(0f, q.x, 0);
        assertEquals(0f, q.y, 0);
        assertEquals(0f, q.z, 0);
    }

    public void testSetIdentity() throws Exception {
        Quaternion4f q = new Quaternion4f(2,3,4,5);
        q.setIdentity();
        assertEquals(1f, q.w, 0);
        assertEquals(0f, q.x, 0);
        assertEquals(0f, q.y, 0);
        assertEquals(0f, q.z, 0);
    }

    public void testFromMatrix_identity() throws Exception {
        Quaternion4f q = new Quaternion4f(2, 3, 4, 5);
        q.fromRotation(
            1, 0, 0,
            0, 1, 0,
            0, 0, 1);
        assertEquals(1, q.w, 0);
        assertEquals(0, q.x, 0);
        assertEquals(0, q.y, 0);
        assertEquals(0, q.z, 0);
    }

    public void testFromMatrix_90x() throws Exception {
        Quaternion4f q = new Quaternion4f(2, 3, 4, 5);
        Matrix3f m = new Matrix3f();
        m.rotateX(m, (float)Math.PI/2);
        q.fromMatrix(m);
        assertEquals((float)Math.sqrt(0.5f), q.w, EPSILON);
        assertEquals((float)Math.sqrt(0.5f), q.x, EPSILON);
        assertEquals(0, q.y, 0);
        assertEquals(0, q.z, 0);
    }

    public void testFromMatrix_90y() throws Exception {
        Quaternion4f q = new Quaternion4f(2, 3, 4, 5);
        Matrix3f m = new Matrix3f();
        m.rotateY(m, (float)Math.PI/2);
        q.fromMatrix(m);
        assertEquals((float)Math.sqrt(0.5f), q.w, EPSILON);
        assertEquals(0, q.x, 0);
        assertEquals((float)Math.sqrt(0.5f), q.y, EPSILON);
        assertEquals(0, q.z, 0);
    }

    public void testFromMatrix_90z() throws Exception {
        Quaternion4f q = new Quaternion4f(2,3,4,5);
        Matrix3f m = new Matrix3f();
        m.rotateZ(m, (float)Math.PI/2);
        q.fromMatrix(m);
        assertEquals((float)Math.sqrt(0.5f), q.w, EPSILON);
        assertEquals(0, q.x, 0);
        assertEquals(0, q.y, 0);
        assertEquals((float)Math.sqrt(0.5f), q.z, EPSILON);
    }

    public void testFromAxisAngle_90x() throws Exception {
        Quaternion4f q = new Quaternion4f(2,3,4,5);
        q.fromAxisAngle(1, 0, 0, (float)((float)Math.PI / 2));
        assertEquals( (float)(float)Math.sqrt(0.5f), q.w, EPSILON);
        assertEquals( (float)(float)Math.sqrt(0.5f), q.x, EPSILON);
        assertEquals( 0f, q.y, EPSILON);
        assertEquals( 0f, q.z, EPSILON);
    }

    public void testFromAxisAngle_90y() throws Exception {
        Quaternion4f q = new Quaternion4f(2,3,4,5);
        q.fromAxisAngle(0, 1, 0, (float)((float)Math.PI / 2));
        assertEquals( (float)(float)Math.sqrt(0.5f), q.w, EPSILON);
        assertEquals( 0f, q.x, EPSILON);
        assertEquals( (float)(float)Math.sqrt(0.5f), q.y, EPSILON);
        assertEquals( 0f, q.z, EPSILON);
    }

    public void testFromAxisAngle_90z() throws Exception {
        Quaternion4f q = new Quaternion4f(2,3,4,5);
        q.fromAxisAngle(0, 0, 1, (float)((float)Math.PI / 2));
        assertEquals( (float)(float)Math.sqrt(0.5f), q.w, EPSILON);
        assertEquals( 0, q.x, EPSILON);
        assertEquals( 0, q.y, EPSILON);
        assertEquals( (float)(float)Math.sqrt(0.5f), q.z, EPSILON);
    }

//    public void testDot() throws Exception {
//
//    }
//
//    public void testMul() throws Exception {
//
//    }
//
//    public void testNegate() throws Exception {
//
//    }
//
//    public void testConj() throws Exception {
//
//    }

    public void testTransform_identity() throws Exception {
        Quaternion4f q = new Quaternion4f();
        Vec3f r = new Vec3f();
        float d = (float)Math.sqrt(1 + 2*2 + 3*3);
        r.transform(q, new Vec3f(1/d,2/d,3/d));
        assertEquals(1/d, r.x, EPSILON);
        assertEquals(2/d, r.y, EPSILON);
        assertEquals(3/d, r.z, EPSILON);
    }

    public void testTransform_90x() throws Exception {
        // rotation around x axis by 90 degrees
        Quaternion4f q = new Quaternion4f((float)(float)Math.sqrt(0.5f), (float)-(float)Math.sqrt(0.5f), 0, 0);
        Vec3f r = new Vec3f();
        r.transform(q, new Vec3f(0,0,1));
        assertEquals(0, r.x, EPSILON);
        assertEquals(1, r.y, EPSILON);
        assertEquals(0, r.z, EPSILON);

        r.transform(q, r);
        assertEquals(0, r.x, EPSILON);
        assertEquals(0, r.y, EPSILON);
        assertEquals(-1, r.z, EPSILON);

        r.transform(q, r);
        assertEquals(0, r.x, EPSILON);
        assertEquals(-1, r.y, EPSILON);
        assertEquals(0, r.z, EPSILON);

        r.transform(q, r);
        assertEquals(0, r.x, EPSILON);
        assertEquals(0, r.y, EPSILON);
        assertEquals(1, r.z, EPSILON);
    }

    public void testTransform_90y() throws Exception {

    }

    public void testArcLength() throws Exception {
        Quaternion4f q = new Quaternion4f(0.0f, 1.0f, 0.0f, 0.0f);
        Quaternion4f p = new Quaternion4f(0.0f, 0.0f, 0.0f, 1.0f);
        assertEquals((float)Math.PI, Quaternion4f.arcLength(p, q), EPSILON);

        q.fromAxisAngle(1.0f, 2.0f, 3.0f, 90 * (float)Math.PI / 180);
        p.fromAxisAngle(1.0f, 2.0f, 3.0f, -20 * (float)Math.PI / 180);
        assertEquals(90 + 20, Quaternion4f.arcLength(p, q) * 180 / (float)Math.PI, EPSILON);
    }

    public void testRandom() throws Exception {
        final int N = 1000;
        Random rand = new Random(1);
        Quaternion4f q = new Quaternion4f();
        Quaternion4f p = new Quaternion4f();
        for (int i=0 ; i<N ; ++i) {
            q.random(rand);
            // the result should be normalized
            assertEquals(1.0f, q.dot(q), EPSILON);
            // the result should be different than the previoud
            assertTrue(Quaternion4f.arcLength(q, p) > EPSILON);
            Quaternion4f t = p; p = q; q = t;
        }
    }

    public void testFromMatrix() {
        Quaternion4f q = new Quaternion4f(2,3,4,5);
        Matrix3f m = new Matrix3f();
        q.fromMatrix(m);
        assertQuaternion(q, 1, 0, 0, 0);

        m.rotate(m, 1, 0, 0, (float)Math.PI/6);
        q.fromMatrix(m);
        assertQuaternion(q,
            (float)Math.cos((float)Math.PI/6/2),
            (float)Math.sin((float)Math.PI/6/2), 0, 0);
        Matrix3f n = new Matrix3f();
        n.fromQuaternion(q);
        Matrix3fTest.assertRotation(n,
            m.m00, m.m01, m.m02,
            m.m10, m.m11, m.m12,
            m.m20, m.m21, m.m22);
    }

    public void testQuaternionToFromMatrix() {
        Quaternion4f q = new Quaternion4f();
        Matrix3f m = new Matrix3f();
        Quaternion4f p = new Quaternion4f();
        Random random = new Random(1);
        final int N = 1000;
        for (int i=0 ; i<N ; ++i) {
            q.random(random);
            m.fromQuaternion(q);
            p.fromMatrix(m);
            assertQuaternion(p, q.w, q.x, q.y, q.z);
        }
    }

    public void testEuler() {
        Quaternion4f q = new Quaternion4f();
        q.fromRPY(0, 0, 0);
        assertQuaternion(
            q,
            1,
            0,
            0,
            0);

        q.fromRPY((float)Math.PI / 6, 0, 0);
        assertQuaternion(
            q,
            (float)Math.cos((float)Math.PI / 12),
            (float)Math.sin((float)Math.PI / 12),
            0,
            0);

        q.fromRPY(0, (float)Math.PI / 6, 0);
        assertQuaternion(
            q,
            (float)Math.cos((float)Math.PI / 12),
            0,
            (float)Math.sin((float)Math.PI / 12),
            0);

        q.fromRPY(0, 0, (float)((float)Math.PI / 6));
        assertQuaternion(
            q,
            (float)(float)Math.cos((float)Math.PI / 12),
            0,
            0,
            (float)(float)Math.sin((float)Math.PI / 12));

//        Quaternion4f a = new Quaternion4f();
//        final int N = 1000;
//        final Random random = new Random(1);
//        for (int i = 0; i < N; ++i) {
//            float roll = (random.nextfloat() - 0.5ff) * (float)(2*(float)Math.PI);
//            float pitch = (random.nextfloat() - 0.5ff) * (float)(2*(float)Math.PI);
//            float yaw = (random.nextfloat() - 0.5ff) * (float)(2*(float)Math.PI);
//
//            q.fromRPY(roll, pitch, yaw);
//            a.fromRPYAlt(roll, pitch, yaw);
//
//            assertQuaternion(q, a);
//        }
    }

//    public void testBenchmark() throws Exception {
//        final int WARMUP_CYCLES = 10;
//        int n = 10000;
//        Quaternion4f q = new float[4];
//        Random random = new Random(1);
//        ReferenceQueue<Quaternion4f> refQueue = new ReferenceQueue<>();
//        for (int cycle=0 ;; ++cycle) {
//            Quaternion4f input = new float[n * 3];
//            random.setSeed(1);
//            for (int i = 0; i < n*3; ++i) {
//                input[i] = 2 * (float)Math.PI * (random.nextfloat() - 0.5f);
//            }
//
//            final long timerStart = System.nanoTime();
//            benchmarkFromRPY(n, q, input);
//            final long timerLap = System.nanoTime();
//            benchmarkFromRPYAlt(n, q, input);
//            final long timerEnd = System.nanoTime();
//
//            if (cycle < WARMUP_CYCLES) {
//                continue;
//            }
//            if (timerEnd - timerStart < 100_000_000) {
//                n = (n * 150_000_000L / (timerEnd - timerStart));
//                System.out.ln("Setting N = "+n);
//            } else {
//                System.out.printf("Quaternion.fromRPY    %.0f ns%n",
//                    (timerLap - timerStart) / (float)n);
//                System.out.printf("Quaternion.fromRPYAlt %.0f ns%n",
//                    (timerEnd - timerLap) / (float)n);
//                break;
//            }
//        }
//    }
//
//    private void benchmarkFromRPYAlt(int n, Quaternion4f q, Quaternion4f input) {
//        for (int i = 0; i< n *3 ; i+=3) {
//            Quaternion.fromRPYAlt(q, input[i], input[i+1], input[i+2]);
//        }
//    }
//
//    private void benchmarkFromRPY(int n, Quaternion4f q, Quaternion4f input) {
//        for (int i = 0; i< n *3 ; i+=3) {
//            Quaternion.fromRPY(q, input[i], input[i+1], input[i+2]);
//        }
//    }

//    public void testTransform() {
//        Random random = new Random(1);
//        Quaternion4f q = new Quaternion4f();
//        float[] m = new float[9];
//        float[] v = new float[3];
//        float[] vq = new float[3];
//        float[] vm = new float[3];
//        final int N = 10000;
//        final int M = 1000;
//        for (int i=0 ; i<N ; ++i) {
//            Quaternion.random(q, random);
//            Matrix3.fromQuaternion(m, q);
//            for (int j=0 ; j<M ; ++j) {
//                for (int k=0 ; k<3 ; ++k)
//                    v[k] = random.nextfloat() - 0.5f;
//                Matrix3.transform(vm, m, v);
//                Quaternion.transform(vq, q, v);
//
//                for (int k=0 ; k<3 ; ++k)
//                    assertEquals(vm[k], vq[k], EPSILON);
//            }
//        }
//    }

//    public void testBenchmarkTransform() {
//        Random random = new Random(1);
//        Quaternion4f q = new float[4];
//        float[] r = new float[3];
//
//        Quaternion.random(q, random);
//        float[] v = random.floats().limit(3).toArray();
//        final int N = 10_000_000;
//
//        for (int run=0 ; run<10 ; ++run) {
//            final long timerStart = System.nanoTime();
//            for (int i = 0; i < N; ++i) {
//                Quaternion.transform(r, q, v);
//            }
//            final long timerEnd = System.nanoTime();
//            System.out.printf(
//                "Quaternion transform: %.2f ns%n",
//                (timerEnd - timerStart) / (float) N);
//        }
//    }


//    public void testMul() throws Exception {
//        float[] a = new float[4];
//        float[] b = new float[4];
//        float[] q = new float[4];
//        Random random = new Random();
//        Quaternion.random(a, random);
//        Quaternion.random(b, random);
//        Quaternion.mul(q, a, b);
//        final float[] m = Matrix3.identity();
//        final float[] n = Matrix3.identity();
//        Matrix3.fromQuaternion(m, a);
//        Matrix3.fromQuaternion(n, b);
//        Matrix3.mul(m, m, n);
//        Quaternion.fromMatrix(b, m);
//
//        for (int i=0 ; i<4 ; ++i) {
//            assertEquals();
//        }
//    }
}