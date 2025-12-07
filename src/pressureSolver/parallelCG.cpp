#include "parallelCG.hpp"
#include <mpi.h>
#include <cmath>
#include <limits>
#include <iostream>
#include <cassert>

ParallelCG::ParallelCG(std::shared_ptr<Discretization> discretization, double epsilon,int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning),
    dx2_(discretization_->dx() * discretization_->dx()),
    dy2_(discretization_->dy() * discretization_->dy()),
    direction_(FieldVariable({partitioning_->nCellsLocal()[0] + 3,
        partitioning_->nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()},
        discretization_->meshWidth())),
    residual_(FieldVariable({partitioning_->nCellsLocal()[0] + 3,
        partitioning_->nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()},
        discretization_->meshWidth())),
    w_(FieldVariable({partitioning_->nCellsLocal()[0] + 3,
        partitioning_->nCellsLocal()[1] + 3}, {-1.5 * discretization_->dx(), -1.5 * discretization_->dy()},
        discretization_->meshWidth()))
    {}

void ParallelCG::solve() {
    const double eps_2 = epsilon_ * epsilon_;
    const int N = partitioning_->getNCellsGlobal(); // total number of real cells in all partitions
    int iter = 0;

    const int pIBegin = discretization_->pIBegin();
    const int pIEnd = discretization_->pIEnd();
    const int pJBegin = discretization_->pJBegin();
    const int pJEnd = discretization_->pJEnd();

    residualOld2_ = 0.0;
    const double diag_precond = 1.0 / (2.0 / dx2_ + 2.0 / dy2_); // diagonal preconditioner (scalar)

    // compute initial residual r = b - A*p and set direction = z = M^{-1} r
    for (int i = pIBegin; i <= pIEnd; i++) {
        for (int j = pJBegin; j <= pJEnd; j++) {
            residual_(i,j) = discretization_->rhs(i,j)
                - ((discretization_->p(i+1,j) - 2.0 * discretization_->p(i,j) + discretization_->p(i-1,j)) / dx2_
                 + (discretization_->p(i,j+1) - 2.0 * discretization_->p(i,j) + discretization_->p(i,j-1)) / dy2_);
            direction_(i,j) = residual_(i,j) * diag_precond; // apply scalar preconditioner: z = M^{-1} r
            residualOld2_ += residual_(i,j) * direction_(i,j); // r^T z (local)
        }
    }

    // global sum of initial r^T z
    MPI_Request request_residual;
    MPI_Iallreduce(MPI_IN_PLACE, &residualOld2_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_residual);
    // need direction halos for first A*direction product later
    communicateAndSetBoundaryValuesForDirection();
    MPI_Wait(&request_residual, MPI_STATUS_IGNORE);

    if (residualOld2_ / N < eps_2) {
        return; // initial guess is good enough
    }

    // main CG iteration loop (fused reduction: reduce dTw and wTw together)
    while (iter < maximumNumberOfIterations_ && residualOld2_ / N > eps_2) {
        iter++;

        // compute w = A * direction and the two local dot-products:
        // local_dots[0] = direction^T * w  (dTw_local)
        // local_dots[1] = w^T * w          (wTw_local)
        double local_dots[2] = {0.0, 0.0};
        for (int i = pIBegin; i <= pIEnd; i++) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                const double dir_ij = direction_(i,j);
                const double w_ij = ((direction_(i+1,j) - 2.0 * dir_ij + direction_(i-1,j)) / dx2_
                                  + (direction_(i,j+1) - 2.0 * dir_ij + direction_(i,j-1)) / dy2_);
                w_(i,j) = w_ij;
                local_dots[0] += dir_ij * w_ij; // dTw_local = z^T w
                local_dots[1] += w_ij * w_ij;   // wTw_local = w^T w
            }
        }

        // reduce both scalars in one non-blocking allreduce (fused)
        double global_dots[2] = {0.0, 0.0};
        MPI_Request request_dots;
        MPI_Iallreduce(local_dots, global_dots, 2, MPI_DOUBLE, MPI_SUM, cartComm_, &request_dots);

        // wait for reduction to finish (we must have global dTw to compute alpha)
        MPI_Wait(&request_dots, MPI_STATUS_IGNORE);

        const double dTw_global = global_dots[0];
        const double wTw_global = global_dots[1];

        // guard against zero dTw (break or handle)
        if (dTw_global == 0.0) {
            // breakdown — cannot continue
            break;
        }

        alpha_ = residualOld2_ / dTw_global;

        // compute global residualNew2 using algebra that uses reduced scalars
        // r_new^T z_new = r^T z - 2*alpha*dTw + alpha^2 * diag_precond * (w^T w)
        residualNew2_ = residualOld2_ - 2.0 * alpha_ * dTw_global + alpha_ * alpha_ * diag_precond * wTw_global;

        // compute beta to update search direction
        double beta = 0.0;
        if (residualOld2_ != 0.0) {
            beta = residualNew2_ / residualOld2_;
        } else {
            beta = 0.0;
        }

        // update solution p, residual r, and direction z in-place (use global alpha and beta)
        for (int i = pIBegin; i <= pIEnd; i++) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                discretization_->p(i,j) += alpha_ * direction_(i,j);
                residual_(i,j) -= alpha_ * w_(i,j);
                const double temp_ij = residual_(i,j) * diag_precond; // z_new = M^{-1} r_new
                direction_(i,j) = temp_ij + beta * direction_(i,j);
            }
        }

        // exchange halos for updated direction before next A*direction
        communicateAndSetBoundaryValuesForDirection();

        // update residualOld2_ for next iteration
        residualOld2_ = residualNew2_;
    }

    // final halo exchange for p (solution)
    communicateAndSetBoundaryValues();
}



void ParallelCG::communicateAndSetBoundaryValuesForDirection() {
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
            direction_(i, pJEnd + 1) = direction_(i, pJEnd);
        }
    } else {
        for (int i = pIBegin; i <= pIEnd; i++) {
            sendBufferTop[i - pIBegin] = direction_(i, pJEnd);
        }
        // instantiate non-blocking sends and receives
        MPI_Isend(sendBufferTop.data(), sendBufferTop.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, cartComm_, &requestTop);
        MPI_Irecv(sendBufferTop.data(), sendBufferTop.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, cartComm_, &requestTop); 
        // override request top with the receive request. and it does not matter since we only wait for receives later (if a receive is not done yet, the send cannot be done either)
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            direction_(i, pJBegin - 1) = direction_(i, pJBegin);
        }
    } else {
        for (int i = pIBegin; i <= pIEnd; i++) {
            sendBufferBottom[i - pIBegin] = direction_(i, pJBegin);
        }
        MPI_Isend(sendBufferBottom.data(), sendBufferBottom.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, cartComm_, &requestBottom);
        MPI_Irecv(sendBufferBottom.data(), sendBufferBottom.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, cartComm_, &requestBottom);
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = pJBegin - 1; j <= pJEnd + 1; j++) {
            direction_(pIBegin - 1, j) = direction_(pIBegin, j);
        }
    } else {
        for (int j = pJBegin; j <= pJEnd; j++) {
            sendBufferLeft[j - pJBegin] = direction_(pIBegin, j);
        }
        MPI_Isend(sendBufferLeft.data(), sendBufferLeft.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, cartComm_, &requestLeft);
        MPI_Irecv(sendBufferLeft.data(), sendBufferLeft.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, cartComm_, &requestLeft);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = pJBegin - 1; j <= pJEnd + 1; j++) {
            direction_(pIEnd + 1, j) = direction_(pIEnd, j);
        }
    } else {
        for (int j = pJBegin; j <= pJEnd; j++) {
            sendBufferRight[j - pJBegin] = direction_(pIEnd, j);
        }
        MPI_Isend(sendBufferRight.data(), sendBufferRight.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, cartComm_, &requestRight);
        MPI_Irecv(sendBufferRight.data(), sendBufferRight.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, cartComm_, &requestRight);
    }

    // nachdem kommuniziert wurde setzren wir die ghost nodes mit den empfangenen werten für ränder die nicht am globalen rand liegen
    
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestTop, MPI_STATUS_IGNORE);  // Wait for receive
        for (int i = pIBegin; i <= pIEnd; i++) {
            // set ghost cells
            direction_(i, pJEnd + 1) = sendBufferTop[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestBottom, MPI_STATUS_IGNORE);  // Wait for receive
        for (int i = pIBegin; i <= pIEnd; i++) {
            // set ghost cells
            direction_(i, pJBegin - 1) = sendBufferBottom[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestLeft, MPI_STATUS_IGNORE);  // Wait for receive
        for (int j = pJBegin; j <= pJEnd; j++) {
            // set ghost cells
            direction_(pIBegin - 1, j) = sendBufferLeft[j - pJBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestRight, MPI_STATUS_IGNORE);  // Wait for receive
        for (int j = pJBegin; j <= pJEnd; j++) {
            // set ghost cells
            direction_(pIEnd + 1, j) = sendBufferRight[j - pJBegin];
        }
    }
}