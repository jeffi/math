package edu.unc.cs.robotics.math;

public interface DoubleMatrix extends Matrix {
    double getCoeff(int r, int c);
    void setCoeff(int r, int c, double value);
}
