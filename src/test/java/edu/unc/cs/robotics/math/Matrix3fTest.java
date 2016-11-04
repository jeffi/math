// Automatically generated from Matrix3dTest.java, DO NOT EDIT!
package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/4/16.
 */
public class Matrix3fTest extends TestCase {

    private static void assertMatrix(Matrix3f actual, Matrix3f expected) {
        for (int r=0 ; r<3 ; ++r) {
            for (int c=0; c<3 ; ++c) {
                assertEquals(r+","+c, expected.getCoeff(r,c), actual.getCoeff(r,c), 1e-6f);
            }
        }
    }

    private static void assertMatrix(Matrix3f actual, float ... expectedCoeffs) {
        Matrix3f expected = new Matrix3f(expectedCoeffs).transpose();
        assertMatrix(actual, expected);
    }

    static void assertRotation(
        Matrix3f actual,
        float ... expectedCoeffs)
    {
        assertMatrix(actual, expectedCoeffs);
        assertTrue(actual.isSpecialOrthogonal());
    }

    public void testTranspose() throws Exception {
        Matrix3f arg = new Matrix3f(
            1, 4, 7,
            2, 5, 8,
            3, 6, 9
        );
        Matrix3f r = new Matrix3f();
        r.transpose(arg);
        arg.transpose(arg);
        for (int i=0 ; i<9 ; ++i) {
            assertEquals(i+1.0f, r.getCoeff(i/3,i%3));
            assertEquals(i+1.0f, arg.getCoeff(i/3,i%3));
        }
    }

    public void testRotation() throws Exception {
        Matrix3f r = new Matrix3f();

        float a = ((float)Math.PI/6);
        r.rotation(1, 0, 0, a);
        assertRotation(r,
            1, 0, 0,
            0, (float)Math.cos(a), -(float)Math.sin(a),
            0, (float)Math.sin(a), (float)Math.cos(a));

        r.rotation(0, 1, 0, a);
        assertRotation(r,
            (float)Math.cos(a), 0, (float)Math.sin(a),
            0, 1, 0,
            -(float)Math.sin(a), 0, (float)Math.cos(a));

        r.rotation(0, 0, 1, a);
        assertRotation(r,
            (float)Math.cos(a), -(float)Math.sin(a), 0,
            (float)Math.sin(a), (float)Math.cos(a), 0,
            0, 0, 1);

        r.rotation(-3, 0, 0, a);
        assertRotation(r,
            1, 0, 0,
            0, (float)Math.cos(a), (float)Math.sin(a),
            0,-(float)Math.sin(a), (float)Math.cos(a));

        r.rotation(1, 1, 1, ((float)Math.PI*2/3));
        assertRotation(r,
            0, 0, 1,
            1, 0, 0,
            0, 1, 0);
    }

    public void testRotateX() throws Exception {
        Matrix3f ident = new Matrix3f();
        Matrix3f r = new Matrix3f();

        for (int i=0 ; i<2 ; ++i) {

            if (i == 0)
                r.rotate(ident, 1, 0, 0, 0);
            else
                r.rotateX(ident, 0);

            assertRotation(r,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1);

            if (i == 0)
                r.rotate(ident, 1, 0, 0, ((float)Math.PI / 2));
            else
                r.rotateX(ident, ((float)Math.PI / 2));

            assertRotation(r,
                1, 0, 0,
                0, 0, -1,
                0, 1, 0);

            if (i == 0)
                r.rotate(ident, 1, 0, 0, (-(float)Math.PI / 2));
            else
                r.rotateX(ident, (-(float)Math.PI / 2));

            assertRotation(r,
                1, 0, 0,
                0, 0, 1,
                0,-1, 0);

            if (i == 0)
                r.rotate(ident, 1, 0, 0, ((float)Math.PI / 6));
            else
                r.rotateX(ident, ((float)Math.PI / 6));

            assertRotation(r,
                1, 0, 0,
                0, (float)Math.sqrt(0.75f), -0.5f,
                0, 0.5f, (float)Math.sqrt(0.75f));
        }
    }

    public void testRotateY() throws Exception {
        Matrix3f ident = new Matrix3f();
        Matrix3f r = new Matrix3f();

        for (int i=0 ; i<2 ; ++i) {
            if (i == 0)
                r.rotate(ident, 0, 1, 0, 0);
            else
                r.rotateY(ident, 0);

            assertRotation(r,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1);

            if (i == 0)
                r.rotate(ident, 0, 1, 0, ((float)Math.PI / 2));
            else
                r.rotateY(ident, ((float)Math.PI / 2));

            assertRotation(r,
                0, 0, 1,
                0, 1, 0,
                -1, 0, 0);

            if (i == 0)
                r.rotate(ident, 0, 1, 0, (-(float)Math.PI / 2));
            else
                r.rotateY(ident, (-(float)Math.PI / 2));

            assertRotation(r,
                0, 0, -1,
                0, 1, 0,
                1, 0, 0);

            if (i == 0)
                r.rotate(ident, 0, 1, 0, ((float)Math.PI / 6));
            else
                r.rotateY(ident, ((float)Math.PI / 6));

            assertRotation(r,
                (float)Math.sqrt(0.75f), 0, 0.5f,
                0, 1, 0,
                -0.5f, 0, (float)Math.sqrt(0.75f));
        }
    }

    public void testRotateZ() throws Exception {
        Matrix3f ident = new Matrix3f();
        Matrix3f r =new Matrix3f();

        for (int i=0 ; i<2 ; ++i) {
            if (i == 0)
                r.rotate(ident, 0, 0, 1, 0);
            else
                r.rotateZ(ident, 0);
            assertRotation(r,
                1, 0, 0,
                0, 1, 0,
                0, 0, 1);

            if (i == 0)
                r.rotate(ident, 0, 0, 1, ((float)Math.PI / 2));
            else
                r.rotateZ(ident, ((float)Math.PI / 2));
            assertRotation(r,
                0, -1, 0,
                1, 0, 0,
                0, 0, 1);

            if (i == 0)
                r.rotate(ident, 0, 0, 1, (-(float)Math.PI / 2));
            else
                r.rotateZ(ident, (-(float)Math.PI/2));
            assertRotation(r,
                0, 1, 0,
                -1, 0, 0,
                0, 0, 1);

            if (i == 0)
                r.rotate(ident, 0, 0, 1, ((float)Math.PI / 6));
            else
                r.rotateZ(ident, ((float)Math.PI/6));
            assertRotation(r,
                (float)Math.sqrt(0.75f),-0.5f, 0,
                0.5f, (float)Math.sqrt(0.75f), 0,
                0, 0, 1);
        }
    }

    public void testMul() throws Exception {
        Matrix3f a = new Matrix3f();
        Matrix3f b = new Matrix3f();
        Matrix3f r = new Matrix3f(2,3,4,5,6,7,8,9,10);

        r.mul(a, b);
        assertMatrix(r,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1);

        a = new Matrix3f(new float[] {1,4,7,2,5,8,3,6,9});
        r.mul(a, a);
        assertMatrix(r,
            30, 36, 42,
            66, 81, 96,
            102, 126, 150);
        b = new Matrix3f(new float[] {8,2,7,9,3,4,1,6,5});
        r.mul(a, b);
        assertMatrix(r,
            33, 27, 28,
            84, 75, 64,
            135, 123, 100);
    }
}