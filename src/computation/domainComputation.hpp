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
#include "pressureSolver/domainCG.hpp"
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
         * it is sufficient to only go though all cells and only look at the right and top face since then we will go through all faces exactly once. 
         * the values set here are only when we have dirichlet BCs directly orthogonally flowing through the face direction. 
         * These are onyl set once in the beginning and then never touched again. 
         * we set also the mirror values for the parralel velocities at the faces here to make the first calculation of the time step witdh correct.
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
        std::vector<CellInfo> fluidGhostCellsInfoList_;
        std::vector<CellInfo> fluidCellsInfoList_;
        std::vector<double> sendBufferTopU_;
        std::vector<double> sendBufferTopV_;
        std::vector<double> sendBufferBottomU_;
        std::vector<double> sendBufferBottomV_;
        std::vector<double> sendBufferLeftU_;
        std::vector<double> sendBufferLeftV_;
        std::vector<double> sendBufferRightU_;
        std::vector<double> sendBufferRightV_;
};
