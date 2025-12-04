#pragma once

#include "computation.hpp"
#include "partitioning/partitioning.hpp"
#include <mpi.h>
#include "output_writer/output_writer_paraview_parallel.hpp"
#include "output_writer/output_writer_text_parallel.hpp"
#include "pressureSolver/parallelPressureSolver.hpp"
#include "pressureSolver/RedBlackGaussSeidel.hpp"
#include "pressureSolver/RedBlackSOR.hpp"

class ParallelComputation : public Computation {
    public:
        /**
         * Initialize the parallel computation with settings from file
         */
        void initialize(int argc, char *argv[]) override;

        /**
         * Run the main simulation loop in parallel
         */
        void runSimulation() override;

    protected:
        /**
         * Apply and Communicate boundary values to the ghost nodes
         */
        void applyBoundaryValues() override;
        
        /**
         * Apply initial boundary values to the ghost nodes before starting the simulation
         */
        void applyInitialBoundaryValues() override;

        /**
         * Compute the time step width dt based on the CFL condition and diffusion limits in parallel
         */
        void computeTimeStepWidth() override;

        std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaview_;
        std::unique_ptr<OutputWriterTextParallel> outputWriterText_;
};
