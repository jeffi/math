package edu.unc.cs.robotics.math;

public interface FloatMatrix extends Matrix {
    float getCoeff(int r, int c);
    void setCoeff(int r, int c, float value);
}
