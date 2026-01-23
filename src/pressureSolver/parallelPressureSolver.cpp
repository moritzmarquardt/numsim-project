#include "pressureSolver/parallelPressureSolver.hpp"
#include "pressureSolver/PressureSolver.hpp"
#include "discretization/discretization.hpp"
#include "partitioning/partitioning.hpp"
#include <cmath>

ParallelPressureSolver::ParallelPressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning) {
    cartComm_ = partitioning_->getCartComm();
}

void ParallelPressureSolver::computeResidualNorm() {
    double residual_norm_squared = 0.0;
    int N_local = 0;
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();

    for (int i = discretization_->pIBegin(); i <= discretization_->pIEnd(); i++) {
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
            if (!discretization_->isActivePressureCell(i, j)) {
                continue;
            }
            N_local++;
            const double rhs_ij = discretization_->rhs(i,j);
            const double p_ij = discretization_->p(i,j);
            const double p_e = discretization_->isFluidCell(i+1,j) ? discretization_->p(i+1,j) : p_ij;
            const double p_w = discretization_->isFluidCell(i-1,j) ? discretization_->p(i-1,j) : p_ij;
            const double p_n = discretization_->isFluidCell(i,j+1) ? discretization_->p(i,j+1) : p_ij;
            const double p_s = discretization_->isFluidCell(i,j-1) ? discretization_->p(i,j-1) : p_ij;
            const double p_xx = (p_e - 2.0 * p_ij + p_w) / dx2;
            const double p_yy = (p_n - 2.0 * p_ij + p_s) / dy2;
            const double residual_ij = rhs_ij - (p_xx + p_yy);
            residual_norm_squared += residual_ij * residual_ij;
        }
    }

    double residual_norm_global = 0.0;
    int N_global = 0;
    // perform global reduction to sum up residual norms from all processes
    // blocking call is acceptable here since we need the result before proceeding
    MPI_Allreduce(&residual_norm_squared, &residual_norm_global, 1, MPI_DOUBLE, MPI_SUM, cartComm_);
    MPI_Allreduce(&N_local, &N_global, 1, MPI_INT, MPI_SUM, cartComm_);
    if (N_global > 0) {
        residualNorm_ = residual_norm_global / N_global;
    } else {
        residualNorm_ = 0.0;
    }
}

void ParallelPressureSolver::communicateAndSetBoundaryValues() {

    //TODO: Idea for improvement: optimize setting of boundary values by only sending/receiving the necessary values instead of the whole rows/columns (only red boundary values needed or black)
    
    const int pIBegin = discretization_->pIBegin();
    const int pIEnd = discretization_->pIEnd();
    const int pJBegin = discretization_->pJBegin();
    const int pJEnd = discretization_->pJEnd();

    // buffers for sending and receiving data
    std::vector<double> sendBufferTop(pIEnd - pIBegin + 1, 0.0);
    std::vector<double> sendBufferBottom(pIEnd - pIBegin + 1, 0.0);
    std::vector<double> sendBufferLeft(pJEnd - pJBegin + 1, 0.0);
    std::vector<double> sendBufferRight(pJEnd - pJBegin + 1, 0.0);

    // init MPI request variables
    MPI_Request requestTop, requestBottom, requestLeft, requestRight;

    // Fill easy boundary conditions that do not need communication
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            discretization_->p(i, pJEnd + 1) = discretization_->p(i, pJEnd);
        }
    } else {
        for (int i = pIBegin; i <= pIEnd; i++) {
            sendBufferTop[i - pIBegin] = discretization_->p(i, pJEnd);
        }
        // instantiate non-blocking sends and receives
        MPI_Isend(sendBufferTop.data(), sendBufferTop.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, cartComm_, &requestTop);
        MPI_Irecv(sendBufferTop.data(), sendBufferTop.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, cartComm_, &requestTop); 
        // override request top with the receive request. and it does not matter since we only wait for receives later (if a receive is not done yet, the send cannot be done either)
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            discretization_->p(i, pJBegin - 1) = discretization_->p(i, pJBegin);
        }
    } else {
        for (int i = pIBegin; i <= pIEnd; i++) {
            sendBufferBottom[i - pIBegin] = discretization_->p(i, pJBegin);
        }
        MPI_Isend(sendBufferBottom.data(), sendBufferBottom.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, cartComm_, &requestBottom);
        MPI_Irecv(sendBufferBottom.data(), sendBufferBottom.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, cartComm_, &requestBottom);
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = pJBegin - 1; j <= pJEnd + 1; j++) {
            discretization_->p(pIBegin - 1, j) = discretization_->p(pIBegin, j);
        }
    } else {
        for (int j = pJBegin; j <= pJEnd; j++) {
            sendBufferLeft[j - pJBegin] = discretization_->p(pIBegin, j);
        }
        MPI_Isend(sendBufferLeft.data(), sendBufferLeft.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, cartComm_, &requestLeft);
        MPI_Irecv(sendBufferLeft.data(), sendBufferLeft.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, cartComm_, &requestLeft);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = pJBegin - 1; j <= pJEnd + 1; j++) {
            discretization_->p(pIEnd + 1, j) = discretization_->p(pIEnd, j);
        }
    } else {
        for (int j = pJBegin; j <= pJEnd; j++) {
            sendBufferRight[j - pJBegin] = discretization_->p(pIEnd, j);
        }
        MPI_Isend(sendBufferRight.data(), sendBufferRight.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, cartComm_, &requestRight);
        MPI_Irecv(sendBufferRight.data(), sendBufferRight.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, cartComm_, &requestRight);
    }

    // nachdem kommuniziert wurde setzren wir die ghost nodes mit den empfangenen werten für ränder die nicht am globalen rand liegen
    
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);  // Wait for receive
        for (int i = pIBegin; i <= pIEnd; i++) {
            // set ghost cells
            discretization_->p(i, pJEnd + 1) = sendBufferTop[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);  // Wait for receive
        for (int i = pIBegin; i <= pIEnd; i++) {
            // set ghost cells
            discretization_->p(i, pJBegin - 1) = sendBufferBottom[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);  // Wait for receive
        for (int j = pJBegin; j <= pJEnd; j++) {
            // set ghost cells
            discretization_->p(pIBegin - 1, j) = sendBufferLeft[j - pJBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);  // Wait for receive
        for (int j = pJBegin; j <= pJEnd; j++) {
            // set ghost cells
            discretization_->p(pIEnd + 1, j) = sendBufferRight[j - pJBegin];
        }
    }
}
