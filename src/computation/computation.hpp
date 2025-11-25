#pragma once
#include "settings.hpp"
#include "discretization/discretization.hpp"
#include "discretization/centralDifferences.hpp"
#include "discretization/donorCell.hpp"
#include "pressureSolver/PressureSolver.hpp"
#include "pressureSolver/GaussSeidel.hpp"
#include "pressureSolver/SOR.hpp"
#include "output_writer/output_writer_paraview.hpp"
#include "output_writer/output_writer_text.hpp"


class Computation {
public:
    virtual void initialize(int argc, char *argv[]);
    virtual void runSimulation();

protected:
    virtual void computeTimeStepWidth();
    virtual void applyBoundaryValues();
    virtual void applyInitialBoundaryValues();
    void computePreliminaryVelocities();
    void computeRightHandSide();
    void computePressure();
    void computeVelocities();

    Settings settings_;
    std::shared_ptr<Discretization> discretization_;
    std::unique_ptr<PressureSolver> pressureSolver_;
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    std::unique_ptr<OutputWriterText> outputWriterText_;
    std::array<double, 2> meshWidth_;
    double dt_;
};
