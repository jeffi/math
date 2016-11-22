package edu.unc.cs.robotics.math;

import java.io.Serializable;

/**
 * This is a very incomplete implementation of complex numbers following the
 * same reuse pattern of this package.  Currently sufficient to help with
 * the implementation of the quartic polynomial solver.
 */
public class Complex2d implements Serializable, Cloneable {
    private static final int HASH_BASE = Complex2d.class.getName().hashCode();

    private static final long serialVersionUID = -3859444351252946830L;

    public double real;
    public double imaginary;

    public Complex2d() {
    }

    public Complex2d(double real) {
        this.real = real;
    }

    public Complex2d(double real, double imaginary) {
        this.real = real;
        this.imaginary = imaginary;
    }

    public static double abs(double real, double imaginary) {
        // TODO: handle NaN & Inf
        if (Math.abs(real) < Math.abs(imaginary)) {
            if (imaginary == 0.0) {
                return Math.abs(real);
            }
            double q = real / imaginary;
            return Math.abs(imaginary) * Math.sqrt(1 + q*q);
        } else {
            if (real == 0.0) {
                return Math.abs(imaginary);
            }
            double q = imaginary / real;
            return Math.abs(real) * Math.sqrt(1 + q*q);
        }
    }

    public double abs() {
        return abs(this.real, this.imaginary);
    }

    /**
     * Returns the complex argument for an imaginary number.  This is the
     * counterclockwise angle from the positive real axis.
     *
     * See also <a href="http://mathworld.wolfram.com/ComplexArgument.html">Complex Argument on Mathworld</a>.
     *
     * @param real the real component of a complex number
     * @param imaginary the imaginary component of a complex number
     * @return the complex argument
     */
    public static double arg(double real, double imaginary) {
        return Math.atan2(imaginary, real);
    }

    /**
     * Returns the complex argument for an imaginary number.  This is the
     * counterclockwise angle from the positive real axis.
     *
     * See also <a href="http://mathworld.wolfram.com/ComplexArgument.html">Complex Argument on Mathworld</a>.
     *
     * @return the complex argument
     */
    public double arg() {
        return Math.atan2(imaginary, real);
    }

    public Complex2d add(double v) {
        this.real += v;
        return this;
    }

    public Complex2d add(double real, double imaginary) {
        this.real += real;
        this.imaginary += imaginary;
        return this;
    }

    public Complex2d add(Complex2d a) {
        return this.add(a.real, a.imaginary);
    }

    public Complex2d sub(double real) {
        this.real -= real;
        return this;
    }

    public Complex2d sub(double real, double imaginary) {
        this.real -= real;
        this.imaginary -= imaginary;
        return this;
    }

    public Complex2d sub(Complex2d a) {
        return this.sub(a.real, a.imaginary);
    }

    public Complex2d mul(double v) {
        // TODO: handle NaN and Inf
        this.real *= v;
        this.imaginary *= v;
        return this;
    }

    public Complex2d mul(double ar, double ai, double br, double bi) {
        this.real      = ar * br - ai * bi;
        this.imaginary = ar * bi + ai * br;
        return this;
    }

    public Complex2d mul(double real, double imaginary) {
        return mul(this.real, this.imaginary, real, imaginary);
    }

    public Complex2d mul(Complex2d a) {
        return mul(a.real, a.imaginary);
    }

    public Complex2d div(double v) {
        // TODO: handle NaN and Inf
        this.real /= v;
        this.imaginary /= v;
        return this;
    }

    /**
     * Computes complex division, and stores result in this.
     *
     * @param numReal real component of numerator
     * @param numImaginary imaginary component of numerator
     * @param denReal real component of denominator
     * @param denImaginary imaginary component of denominator
     * @return {@code this}
     */
    public Complex2d div(
        double numReal, double numImaginary,
        double denReal, double denImaginary)
    {
        if (Math.abs(denReal) < Math.abs(denImaginary)) {
            final double q = denReal/denImaginary;
            final double den = denReal*q + denImaginary;
            this.real = (numReal*q + numImaginary)/den;
            this.imaginary = (numImaginary*q - numReal)/den;
        } else {
            final double q = denImaginary/denReal;
            final double den = denImaginary*q + denReal;
            this.real = (numImaginary*q + numReal)/den;
            this.imaginary = (numImaginary - numReal*q)/den;
        }

        return this;
    }

    public Complex2d div(Complex2d a) {
        return div(this.real, this.imaginary, a.real, a.imaginary);
    }

    public Complex2d sqrt(double real, double imaginary) {
        // TODO: handle NaN
        if (real == 0.0 && imaginary == 0.0) {
            this.real = 0.0;
            this.imaginary = 0.0;
        } else {
            final double t = Math.sqrt((Math.abs(real) + abs(real, imaginary))/2.0);
            if (real >= 0.0) {
                this.real = t;
                this.imaginary = imaginary/(2.0*t);
            } else {
                this.real = Math.abs(imaginary)/(2.0*t);
                this.imaginary = Math.copySign(1.0, imaginary)*t;
            }
        }
        return this;
    }

    public Complex2d sqrt(Complex2d a) {
        return sqrt(a.real, a.imaginary);
    }

    public Complex2d sqrt() {
        return sqrt(this);
    }

    public Complex2d[] cbrt() {
        // the next line should probably be a Math.cbrt(),
        // but pow seems faster and is oddly more accurate in my test
        // case.  And we don't need to handle cbrt of a negative.
        final double absCbrt = Math.pow(abs(), 1.0/3.0);
        final double nthPhi = arg()/3.0;
        final double slice = 2.0 * Math.PI / 3.0;
        final Complex2d[] result = new Complex2d[3];
        for (int i=0 ; i<3 ; ++i) {
            result[i] = new Complex2d(
                absCbrt*Math.cos(nthPhi + i*slice),
                absCbrt*Math.sin(nthPhi + i*slice));
        }
        return result;
    }

    @Override
    public Complex2d clone() {
        try {
            return (Complex2d)super.clone();
        } catch (CloneNotSupportedException e) {
            throw new AssertionError(e);
        }
    }

    @Override
    public int hashCode() {
        long h = Double.doubleToLongBits(real);
        h = h*31 + Double.doubleToLongBits(imaginary);
        return (int)h ^ (int)(h >>> 32) ^ HASH_BASE;
    }

    @Override
    public boolean equals(Object obj) {
        if (obj == this) return true;
        if (obj == null || obj.getClass() != getClass()) return false;
        final Complex2d that = (Complex2d)obj;
        return MathEx.bitEquals(this.real, that.real)
               && MathEx.bitEquals(this.imaginary, that.imaginary);
    }

    @Override
    public String toString() {
        return real + (imaginary < 0 ? imaginary+"i" : "+"+imaginary+"i");
    }
}
