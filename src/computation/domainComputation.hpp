#pragma once

#include "computation.hpp"
#include "partitioning/partitioning.hpp"
#include <mpi.h>
#include "output_writer/output_writer_paraview_parallel.hpp"
#include "output_writer/output_writer_text_parallel.hpp"
#include "pressureSolver/parallelPressureSolver.hpp"
#include "pressureSolver/RedBlackGaussSeidel.hpp"
#include "pressureSolver/domainPressureSolver.hpp"
#include "pressureSolver/domainRBGaussSeidel.hpp"
#include "pressureSolver/RedBlackSOR.hpp"
#include "pressureSolver/parallelCG.hpp"
#include "domain/domain.hpp"

class DomainComputation : public Computation {
    public:
        /**
         * Initialize the parallel computation with settings from file
         */
        void initialize(int argc, char *argv[]) override;

        /**
         * Run the main simulation loop in parallel
         */
        void runSimulation() override;

        void printProgress(double &currentTime, int &iterationCount);

        /**
         * Get the own MPI rank number
         */
        int getRankNo() const {
            return partitioning_->ownRankNo();
        }

    protected:        
        /**
         * Apply initial boundary values to the ghost nodes before starting the simulation
         */
        void applyInitialBoundaryValues() override;

        void communicateGhostCells();

        void computePreliminaryVelocities() override;

        void computeRightHandSide() override;

        void computePressure() override;
        
        void computeVelocities() override;

        void computeTimeStepWidth() override;

        // all the stencils of dicretisation and donor cell need to be implemented using ghost stencils
        double computeD2uDx2(double u_ip1_j, double u_i_j, double u_im1_j) const;

        double computeD2uDy2(double u_i_jp1, double u_i_j, double u_i_jm1) const;

        double computeD2vDx2(double v_ip1_j, double v_i_j, double v_im1_j) const;

        double computeD2vDy2(double v_i_jp1, double v_i_j, double v_i_jm1) const;

        double computeDpDx(double p_ip1_j, double p_i_j) const;

        double computeDpDy(double p_i_jp1, double p_i_j) const;

        double computeDu2Dx(double u_i_j, double u_im1_j, double u_ip1_j) const;

        double computeDv2Dy(double v_i_j, double v_i_jp1, double v_i_jm1) const;

        double computeDuvDx(double u_i_j, double u_i_jp1, double u_im1_j, double u_im1_jp1,
                            double v_i_j, double v_ip1_j, double v_im1_j) const;

        double computeDuvDy(double u_i_j, double u_i_jp1, double u_i_jm1,
                            double v_i_j, double v_ip1_j, double v_i_jm1, double v_ip1_jm1) const;

        std::unique_ptr<OutputWriterParaviewParallel> outputWriterParaview_;
        std::unique_ptr<OutputWriterTextParallel> outputWriterText_;
        MPI_Comm cartComm_;
        std::shared_ptr<Domain> domain_;
};
