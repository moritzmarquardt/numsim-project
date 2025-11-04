#pragma once

class Computation {
public:
    void initialize(int argc, char *argv[]);
    void runSimulation();

private:
    /**
     * C
     */
    void computeTimeStepWidth();
    void applyBoundaryValues();
    void computePreliminaryVelocities();
    void computeRightHandSide();
    void computePressure();
    void computeVelocities();
};
