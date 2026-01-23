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
#include "partitioning/partitioning.hpp"


class Computation {
public:
    /**
     * Initialize the computation with settings from file
     */
    virtual void initialize(int argc, char *argv[]);

    /**
     * Run the main simulation loop
     */
    virtual void runSimulation();

protected:
    /**
     * Compute the time step width dt based on the CFL condition and diffusion limits
     */
    virtual void computeTimeStepWidth();

    /**
     * Apply boundary values to the ghost nodes at each time step
     */
    virtual void applyBoundaryValues();

    /**
     * Apply initial boundary values to the ghost nodes before starting the simulation
     */
    virtual void applyInitialBoundaryValues();

    /**
     * Compute preliminary velocities F and G
     */
    void computePreliminaryVelocities();

    /**
     * Compute the right-hand side of the pressure Poisson equation
     */
    void computeRightHandSide();

    /**
     * Compute the pressure field by solving the Poisson equation
     */
    void computePressure();

    /**
     * Compute the final velocities u and v
     */
    void computeVelocities();

    Settings settings_;
    int simNumber_;
    std::shared_ptr<Discretization> discretization_;
    std::unique_ptr<PressureSolver> pressureSolver_;
    std::unique_ptr<OutputWriterParaview> outputWriterParaview_;
    std::unique_ptr<OutputWriterText> outputWriterText_;
    std::array<double, 2> meshWidth_;
    double dt_;
    std::shared_ptr<Partitioning> partitioning_;
};
