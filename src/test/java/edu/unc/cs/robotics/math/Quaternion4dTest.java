package edu.unc.cs.robotics.math;

import java.util.Random;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/4/16.
 */
public class Quaternion4dTest extends TestCase {
    static final double MAX_ERROR = 1e-9f;

    private static void assertQuaternion(Quaternion4d actual, double ... expected) {
        boolean posEq = true;
        boolean negEq = true;
        double sum = 0;
        for (int i=0 ; i<4 ; ++i) {
            double v = actual.getCoeff(i);
            posEq &= Math.abs(expected[i] - v) < MAX_ERROR;
            negEq &= Math.abs(expected[i] + v) < MAX_ERROR;
            sum += v*v;
        }
        if (!posEq && !negEq) {
            failNotEquals("", new Quaternion4d(expected), actual);
        }
        assertEquals(1.0, sum, 1e-9);
    }

    public void testCreateIdentity() throws Exception {
        final Quaternion4d q = new Quaternion4d();
        assertEquals(1f, q.w, 0);
        assertEquals(0f, q.x, 0);
        assertEquals(0f, q.y, 0);
        assertEquals(0f, q.z, 0);
    }

    public void testSetIdentity() throws Exception {
        Quaternion4d q = new Quaternion4d(2,3,4,5);
        q.setIdentity();
        assertEquals(1f, q.w, 0);
        assertEquals(0f, q.x, 0);
        assertEquals(0f, q.y, 0);
        assertEquals(0f, q.z, 0);
    }

    public void testFromMatrix_identity() throws Exception {
        Quaternion4d q = new Quaternion4d(2, 3, 4, 5);
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
        Quaternion4d q = new Quaternion4d(2, 3, 4, 5);
        Matrix3d m = new Matrix3d();
        m.rotateX(m, Math.PI/2);
        q.fromMatrix(m);
        assertEquals(Math.sqrt(0.5), q.w, 1e-9);
        assertEquals(Math.sqrt(0.5), q.x, 1e-9);
        assertEquals(0, q.y, 0);
        assertEquals(0, q.z, 0);
    }

    public void testFromMatrix_90y() throws Exception {
        Quaternion4d q = new Quaternion4d(2, 3, 4, 5);
        Matrix3d m = new Matrix3d();
        m.rotateY(m, Math.PI/2);
        q.fromMatrix(m);
        assertEquals(Math.sqrt(0.5), q.w, 1e-9);
        assertEquals(0, q.x, 0);
        assertEquals(Math.sqrt(0.5), q.y, 1e-9);
        assertEquals(0, q.z, 0);
    }

    public void testFromMatrix_90z() throws Exception {
        Quaternion4d q = new Quaternion4d(2,3,4,5);
        Matrix3d m = new Matrix3d();
        m.rotateZ(m, Math.PI/2);
        q.fromMatrix(m);
        assertEquals(Math.sqrt(0.5), q.w, 1e-9);
        assertEquals(0, q.x, 0);
        assertEquals(0, q.y, 0);
        assertEquals(Math.sqrt(0.5), q.z, 1e-9);
    }

    public void testFromAxisAngle_90x() throws Exception {
        Quaternion4d q = new Quaternion4d(2,3,4,5);
        q.fromAxisAngle(1, 0, 0, (double)(Math.PI / 2));
        assertEquals( (double)Math.sqrt(0.5), q.w, MAX_ERROR);
        assertEquals( (double)Math.sqrt(0.5), q.x, MAX_ERROR);
        assertEquals( 0f, q.y, MAX_ERROR);
        assertEquals( 0f, q.z, MAX_ERROR);
    }

    public void testFromAxisAngle_90y() throws Exception {
        Quaternion4d q = new Quaternion4d(2,3,4,5);
        q.fromAxisAngle(0, 1, 0, (double)(Math.PI / 2));
        assertEquals( (double)Math.sqrt(0.5), q.w, MAX_ERROR);
        assertEquals( 0f, q.x, MAX_ERROR);
        assertEquals( (double)Math.sqrt(0.5), q.y, MAX_ERROR);
        assertEquals( 0f, q.z, MAX_ERROR);
    }

    public void testFromAxisAngle_90z() throws Exception {
        Quaternion4d q = new Quaternion4d(2,3,4,5);
        q.fromAxisAngle(0, 0, 1, (double)(Math.PI / 2));
        assertEquals( (double)Math.sqrt(0.5), q.w, MAX_ERROR);
        assertEquals( 0, q.x, MAX_ERROR);
        assertEquals( 0, q.y, MAX_ERROR);
        assertEquals( (double)Math.sqrt(0.5), q.z, MAX_ERROR);
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
        Quaternion4d q = new Quaternion4d();
        Vec3d r = new Vec3d();
        double d = Math.sqrt(1 + 2*2 + 3*3);
        r.transform(q, new Vec3d(1/d,2/d,3/d));
        assertEquals(1/d, r.x, MAX_ERROR);
        assertEquals(2/d, r.y, MAX_ERROR);
        assertEquals(3/d, r.z, MAX_ERROR);
    }

    public void testTransform_90x() throws Exception {
        // rotation around x axis by 90 degrees
        Quaternion4d q = new Quaternion4d((double)Math.sqrt(0.5), (double)-Math.sqrt(0.5), 0, 0);
        Vec3d r = new Vec3d();
        r.transform(q, new Vec3d(0,0,1));
        assertEquals(0, r.x, MAX_ERROR);
        assertEquals(1, r.y, MAX_ERROR);
        assertEquals(0, r.z, MAX_ERROR);

        r.transform(q, r);
        assertEquals(0, r.x, MAX_ERROR);
        assertEquals(0, r.y, MAX_ERROR);
        assertEquals(-1, r.z, MAX_ERROR);

        r.transform(q, r);
        assertEquals(0, r.x, MAX_ERROR);
        assertEquals(-1, r.y, MAX_ERROR);
        assertEquals(0, r.z, MAX_ERROR);

        r.transform(q, r);
        assertEquals(0, r.x, MAX_ERROR);
        assertEquals(0, r.y, MAX_ERROR);
        assertEquals(1, r.z, MAX_ERROR);
    }

    public void testTransform_90y() throws Exception {

    }

    public void testRandom() throws Exception {
        final int N = 1000;
        Random rand = new Random(1);
        Quaternion4d q = new Quaternion4d();
        Quaternion4d p = new Quaternion4d();
        for (int i=0 ; i<N ; ++i) {
            q.random(rand);
            // the result should be normalized
            assertEquals(1.0, q.dot(q), MAX_ERROR);
            // the result should be different than the previoud
            assertTrue(Quaternion4d.arcLength(q, p) > MAX_ERROR);
            Quaternion4d t = p; p = q; q = t;
        }
    }

    public void testFromMatrix() {
        Quaternion4d q = new Quaternion4d(2,3,4,5);
        Matrix3d m = new Matrix3d();
        q.fromMatrix(m);
        assertQuaternion(q, 1, 0, 0, 0);

        m.rotate(m, 1, 0, 0, Math.PI/6);
        q.fromMatrix(m);
        assertQuaternion(q,
            Math.cos(Math.PI/6/2),
            Math.sin(Math.PI/6/2), 0, 0);
        Matrix3d n = new Matrix3d();
        n.fromQuaternion(q);
        Matrix3dTest.assertRotation(n,
            m.m00, m.m01, m.m02,
            m.m10, m.m11, m.m12,
            m.m20, m.m21, m.m22);
    }

    public void testQuaternionToFromMatrix() {
        Quaternion4d q = new Quaternion4d();
        Matrix3d m = new Matrix3d();
        Quaternion4d p = new Quaternion4d();
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
        Quaternion4d q = new Quaternion4d();
        q.fromRPY(0, 0, 0);
        assertQuaternion(
            q,
            1,
            0,
            0,
            0);

        q.fromRPY(Math.PI / 6, 0, 0);
        assertQuaternion(
            q,
            Math.cos(Math.PI / 12),
            Math.sin(Math.PI / 12),
            0,
            0);

        q.fromRPY(0, Math.PI / 6, 0);
        assertQuaternion(
            q,
            Math.cos(Math.PI / 12),
            0,
            Math.sin(Math.PI / 12),
            0);

        q.fromRPY(0, 0, (double)(Math.PI / 6));
        assertQuaternion(
            q,
            (double)Math.cos(Math.PI / 12),
            0,
            0,
            (double)Math.sin(Math.PI / 12));

//        Quaternion4d a = new Quaternion4d();
//        final int N = 1000;
//        final Random random = new Random(1);
//        for (int i = 0; i < N; ++i) {
//            double roll = (random.nextdouble() - 0.5f) * (double)(2*Math.PI);
//            double pitch = (random.nextdouble() - 0.5f) * (double)(2*Math.PI);
//            double yaw = (random.nextdouble() - 0.5f) * (double)(2*Math.PI);
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
//        Quaternion4d q = new double[4];
//        Random random = new Random(1);
//        ReferenceQueue<Quaternion4d> refQueue = new ReferenceQueue<>();
//        for (int cycle=0 ;; ++cycle) {
//            Quaternion4d input = new double[n * 3];
//            random.setSeed(1);
//            for (int i = 0; i < n*3; ++i) {
//                input[i] = 2 * Math.PI * (random.nextdouble() - 0.5);
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
//                n = (int)(n * 150_000_000L / (timerEnd - timerStart));
//                System.out.ln("Setting N = "+n);
//            } else {
//                System.out.printf("Quaternion.fromRPY    %.0f ns%n",
//                    (timerLap - timerStart) / (double)n);
//                System.out.printf("Quaternion.fromRPYAlt %.0f ns%n",
//                    (timerEnd - timerLap) / (double)n);
//                break;
//            }
//        }
//    }
//
//    private void benchmarkFromRPYAlt(int n, Quaternion4d q, Quaternion4d input) {
//        for (int i = 0; i< n *3 ; i+=3) {
//            Quaternion.fromRPYAlt(q, input[i], input[i+1], input[i+2]);
//        }
//    }
//
//    private void benchmarkFromRPY(int n, Quaternion4d q, Quaternion4d input) {
//        for (int i = 0; i< n *3 ; i+=3) {
//            Quaternion.fromRPY(q, input[i], input[i+1], input[i+2]);
//        }
//    }

//    public void testTransform() {
//        Random random = new Random(1);
//        Quaternion4d q = new Quaternion4d();
//        double[] m = new double[9];
//        double[] v = new double[3];
//        double[] vq = new double[3];
//        double[] vm = new double[3];
//        final int N = 10000;
//        final int M = 1000;
//        for (int i=0 ; i<N ; ++i) {
//            Quaternion.random(q, random);
//            Matrix3.fromQuaternion(m, q);
//            for (int j=0 ; j<M ; ++j) {
//                for (int k=0 ; k<3 ; ++k)
//                    v[k] = random.nextdouble() - 0.5;
//                Matrix3.transform(vm, m, v);
//                Quaternion.transform(vq, q, v);
//
//                for (int k=0 ; k<3 ; ++k)
//                    assertEquals(vm[k], vq[k], 1e-9);
//            }
//        }
//    }

//    public void testBenchmarkTransform() {
//        Random random = new Random(1);
//        Quaternion4d q = new double[4];
//        double[] r = new double[3];
//
//        Quaternion.random(q, random);
//        double[] v = random.doubles().limit(3).toArray();
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
//                (timerEnd - timerStart) / (double) N);
//        }
//    }


//    public void testMul() throws Exception {
//        double[] a = new double[4];
//        double[] b = new double[4];
//        double[] q = new double[4];
//        Random random = new Random();
//        Quaternion.random(a, random);
//        Quaternion.random(b, random);
//        Quaternion.mul(q, a, b);
//        final double[] m = Matrix3.identity();
//        final double[] n = Matrix3.identity();
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