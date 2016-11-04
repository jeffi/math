package edu.unc.cs.robotics.math;

import junit.framework.TestCase;

/**
 * Created by jeffi on 11/4/16.
 */
public class Matrix3dTest extends TestCase {

    private static void assertMatrix(Matrix3d actual, Matrix3d expected) {
        for (int r=0 ; r<3 ; ++r) {
            for (int c=0; c<3 ; ++c) {
                assertEquals(r+","+c, expected.getCoeff(r,c), actual.getCoeff(r,c), 1e-6f);
            }
        }
    }

    private static void assertMatrix(Matrix3d actual, double ... expectedCoeffs) {
        Matrix3d expected = new Matrix3d(expectedCoeffs).transpose();
        assertMatrix(actual, expected);
    }

    static void assertRotation(
        Matrix3d actual,
        double ... expectedCoeffs)
    {
        assertMatrix(actual, expectedCoeffs);
        assertTrue(actual.isSpecialOrthogonal());
    }

    public void testTranspose() throws Exception {
        Matrix3d arg = new Matrix3d(
            1, 4, 7,
            2, 5, 8,
            3, 6, 9
        );
        Matrix3d r = new Matrix3d();
        r.transpose(arg);
        arg.transpose(arg);
        for (int i=0 ; i<9 ; ++i) {
            assertEquals(i+1.0, r.getCoeff(i/3,i%3));
            assertEquals(i+1.0, arg.getCoeff(i/3,i%3));
        }
    }

    public void testRotation() throws Exception {
        Matrix3d r = new Matrix3d();

        double a = (Math.PI/6);
        r.rotation(1, 0, 0, a);
        assertRotation(r,
            1, 0, 0,
            0, Math.cos(a), -Math.sin(a),
            0, Math.sin(a), Math.cos(a));

        r.rotation(0, 1, 0, a);
        assertRotation(r,
            Math.cos(a), 0, Math.sin(a),
            0, 1, 0,
            -Math.sin(a), 0, Math.cos(a));

        r.rotation(0, 0, 1, a);
        assertRotation(r,
            Math.cos(a), -Math.sin(a), 0,
            Math.sin(a), Math.cos(a), 0,
            0, 0, 1);

        r.rotation(-3, 0, 0, a);
        assertRotation(r,
            1, 0, 0,
            0, Math.cos(a), Math.sin(a),
            0,-Math.sin(a), Math.cos(a));

        r.rotation(1, 1, 1, (Math.PI*2/3));
        assertRotation(r,
            0, 0, 1,
            1, 0, 0,
            0, 1, 0);
    }

    public void testRotateX() throws Exception {
        Matrix3d ident = new Matrix3d();
        Matrix3d r = new Matrix3d();

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
                r.rotate(ident, 1, 0, 0, (Math.PI / 2));
            else
                r.rotateX(ident, (Math.PI / 2));

            assertRotation(r,
                1, 0, 0,
                0, 0, -1,
                0, 1, 0);

            if (i == 0)
                r.rotate(ident, 1, 0, 0, (-Math.PI / 2));
            else
                r.rotateX(ident, (-Math.PI / 2));

            assertRotation(r,
                1, 0, 0,
                0, 0, 1,
                0,-1, 0);

            if (i == 0)
                r.rotate(ident, 1, 0, 0, (Math.PI / 6));
            else
                r.rotateX(ident, (Math.PI / 6));

            assertRotation(r,
                1, 0, 0,
                0, Math.sqrt(0.75), -0.5f,
                0, 0.5f, Math.sqrt(0.75));
        }
    }

    public void testRotateY() throws Exception {
        Matrix3d ident = new Matrix3d();
        Matrix3d r = new Matrix3d();

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
                r.rotate(ident, 0, 1, 0, (Math.PI / 2));
            else
                r.rotateY(ident, (Math.PI / 2));

            assertRotation(r,
                0, 0, 1,
                0, 1, 0,
                -1, 0, 0);

            if (i == 0)
                r.rotate(ident, 0, 1, 0, (-Math.PI / 2));
            else
                r.rotateY(ident, (-Math.PI / 2));

            assertRotation(r,
                0, 0, -1,
                0, 1, 0,
                1, 0, 0);

            if (i == 0)
                r.rotate(ident, 0, 1, 0, (Math.PI / 6));
            else
                r.rotateY(ident, (Math.PI / 6));

            assertRotation(r,
                Math.sqrt(0.75), 0, 0.5f,
                0, 1, 0,
                -0.5f, 0, Math.sqrt(0.75));
        }
    }

    public void testRotateZ() throws Exception {
        Matrix3d ident = new Matrix3d();
        Matrix3d r =new Matrix3d();

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
                r.rotate(ident, 0, 0, 1, (Math.PI / 2));
            else
                r.rotateZ(ident, (Math.PI / 2));
            assertRotation(r,
                0, -1, 0,
                1, 0, 0,
                0, 0, 1);

            if (i == 0)
                r.rotate(ident, 0, 0, 1, (-Math.PI / 2));
            else
                r.rotateZ(ident, (-Math.PI/2));
            assertRotation(r,
                0, 1, 0,
                -1, 0, 0,
                0, 0, 1);

            if (i == 0)
                r.rotate(ident, 0, 0, 1, (Math.PI / 6));
            else
                r.rotateZ(ident, (Math.PI/6));
            assertRotation(r,
                Math.sqrt(0.75),-0.5f, 0,
                0.5f, Math.sqrt(0.75), 0,
                0, 0, 1);
        }
    }

    public void testMul() throws Exception {
        Matrix3d a = new Matrix3d();
        Matrix3d b = new Matrix3d();
        Matrix3d r = new Matrix3d(2,3,4,5,6,7,8,9,10);

        r.mul(a, b);
        assertMatrix(r,
            1, 0, 0,
            0, 1, 0,
            0, 0, 1);

        a = new Matrix3d(new double[] {1,4,7,2,5,8,3,6,9});
        r.mul(a, a);
        assertMatrix(r,
            30, 36, 42,
            66, 81, 96,
            102, 126, 150);
        b = new Matrix3d(new double[] {8,2,7,9,3,4,1,6,5});
        r.mul(a, b);
        assertMatrix(r,
            33, 27, 28,
            84, 75, 64,
            135, 123, 100);
    }
}