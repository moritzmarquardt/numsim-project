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

// Patched ParallelCG::solve() with safety and correctness fixes
// Requires: <cassert>, <cmath>, <limits>, <iostream>, <mpi.h>
void ParallelCG::solve() {
    const double eps_2 = epsilon_ * epsilon_;
    const int N = partitioning_->getNCellsGlobal(); // total number of real cells in all partitions
    int iter = 0;

    const int pIBegin = discretization_->pIBegin();
    const int pIEnd   = discretization_->pIEnd();
    const int pJBegin = discretization_->pJBegin();
    const int pJEnd   = discretization_->pJEnd();

    residualOld2_ = 0.0;
    const double diag_precond = 1.0 / (2.0 / dx2_ + 2.0 / dy2_); // diagonal used inside GS updates

    // local sizes (interior only)
    const int nI = pIEnd - pIBegin + 1;
    const int nJ = pJEnd - pJBegin + 1;
    const int localSize = nI * nJ;

    // Temporary buffer to keep old direction (p_k) while we compute z_{k+1}
    std::vector<double> dir_old;
    dir_old.resize(localSize);

    int world_rank = 0;
    MPI_Comm_rank(cartComm_, &world_rank);

    // --- preconditioner: block symmetric red-black Gauss-Seidel ---
    auto applyBlockSGSPreconditioner = [&](int num_sweeps) {
        // initial guess: Jacobi
        for (int i = pIBegin; i <= pIEnd; ++i) {
            for (int j = pJBegin; j <= pJEnd; ++j) {
                direction_(i,j) = residual_(i,j) * diag_precond;
            }
        }
        // do sweeps
        for (int sweep = 0; sweep < num_sweeps; ++sweep) {
            // update red (color 0)
            communicateAndSetBoundaryValuesForDirection();
            for (int i = pIBegin; i <= pIEnd; ++i) {
                for (int j = pJBegin; j <= pJEnd; ++j) {
                    if (((i + j) & 1) == 0) { // red
                        const double sum_neigh =
                            (direction_(i+1,j) + direction_(i-1,j)) / dx2_
                          + (direction_(i,j+1) + direction_(i,j-1)) / dy2_;
                        // z_i = (sum_neigh - r_i) / (2/dx2 + 2/dy2)
                        direction_(i,j) = diag_precond * (sum_neigh - residual_(i,j));
                    }
                }
            }
            // update black (color 1)
            communicateAndSetBoundaryValuesForDirection();
            for (int i = pIBegin; i <= pIEnd; ++i) {
                for (int j = pJBegin; j <= pJEnd; ++j) {
                    if (((i + j) & 1) == 1) { // black
                        const double sum_neigh =
                            (direction_(i+1,j) + direction_(i-1,j)) / dx2_
                          + (direction_(i,j+1) + direction_(i,j-1)) / dy2_;
                        direction_(i,j) = diag_precond * (sum_neigh - residual_(i,j));
                    }
                }
            }
        }
        // final halos ready for next A*direction (optional but consistent)
        communicateAndSetBoundaryValuesForDirection();
    };
    // --- end preconditioner ---

    // Ensure solution p has valid halos before computing residual r = b - A*p
    communicateAndSetBoundaryValues();

    // compute initial residual r = b - A*p
    for (int i = pIBegin; i <= pIEnd; ++i) {
        for (int j = pJBegin; j <= pJEnd; ++j) {
            residual_(i,j) = discretization_->rhs(i,j)
                - ((discretization_->p(i+1,j) - 2.0 * discretization_->p(i,j) + discretization_->p(i-1,j)) / dx2_
                 + (discretization_->p(i,j+1) - 2.0 * discretization_->p(i,j) + discretization_->p(i,j-1)) / dy2_);
        }
    }

    // apply block SGS preconditioner to form initial direction (z_0)
    const int num_sgs_sweeps = 2; // tunable: 1..4 often used; 2 is a reasonable default
    applyBlockSGSPreconditioner(num_sgs_sweeps);

    // compute initial residualOld2_ = r^T z (local accumulation)
    residualOld2_ = 0.0;
    for (int i = pIBegin; i <= pIEnd; ++i) {
        for (int j = pJBegin; j <= pJEnd; ++j) {
            residualOld2_ += residual_(i,j) * direction_(i,j);
        }
    }

    // global sum (wait before any other communicator calls to avoid ordering issues)
    MPI_Request request_residual;
    MPI_Iallreduce(MPI_IN_PLACE, &residualOld2_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_residual);
    MPI_Wait(&request_residual, MPI_STATUS_IGNORE);
    // now make sure direction halos are set for later A*direction
    communicateAndSetBoundaryValuesForDirection();

    // basic sanity checks
    if (!std::isfinite(residualOld2_)) {
        if (world_rank == 0) std::cerr << "ParallelCG::solve - NaN/Inf in initial residualOld2_: " << residualOld2_ << "\n";
        return;
    }

    if (residualOld2_ / N < eps_2) {
        return; // initial guess is good enough (preconditioned inner product test retained)
    }

    // main CG iteration loop (we keep the two-reduction approach but with stronger preconditioner)
    double dTw_local = 0.0;
    while (iter < maximumNumberOfIterations_ && residualOld2_ / N > eps_2) {
        iter++;

        // save old direction (p_k) into dir_old before it gets overwritten
        int idx = 0;
        for (int i = pIBegin; i <= pIEnd; ++i) {
            for (int j = pJBegin; j <= pJEnd; ++j) {
                dir_old[idx++] = direction_(i,j);
            }
        }
        // safety: ensure we filled exactly localSize
        assert(idx == localSize && "dir_old size mismatch");

        // compute w = A * direction
        dTw_local = 0.0;
        for (int i = pIBegin; i <= pIEnd; ++i) {
            for (int j = pJBegin; j <= pJEnd; ++j) {
                const double dir_ij = direction_(i,j);
                const double w_ij = ((direction_(i+1,j) - 2.0 * dir_ij + direction_(i-1,j)) / dx2_
                                   + (direction_(i,j+1) - 2.0 * dir_ij + direction_(i,j-1)) / dy2_);
                w_(i,j) = w_ij;
                dTw_local += dir_ij * w_ij;
            }
        }

        // global reduction for dTw (wait before other comms)
        double dTw_global = 0.0;
        MPI_Request request_dTw;
        MPI_Iallreduce(&dTw_local, &dTw_global, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_dTw);
        MPI_Wait(&request_dTw, MPI_STATUS_IGNORE);

        // robust check for breakdown or bad values
        const double denom_tol = std::numeric_limits<double>::min() * 1e6; // tiny absolute threshold
        if (!std::isfinite(dTw_global) || dTw_global <= denom_tol) {
            if (world_rank == 0) {
                std::cerr << "ParallelCG::solve - breakdown: d^T A d = " << dTw_global << " (tol=" << denom_tol << ")\n";
            }
            break;
        }

        // further safety: residualOld2_ must be finite
        if (!std::isfinite(residualOld2_)) {
            if (world_rank == 0) std::cerr << "ParallelCG::solve - nonfinite residualOld2_ before alpha computation\n";
            break;
        }

        alpha_ = residualOld2_ / dTw_global;
        if (!std::isfinite(alpha_)) {
            if (world_rank == 0) std::cerr << "ParallelCG::solve - nonfinite alpha: " << alpha_ << "\n";
            break;
        }

        // update solution p and residual r
        for (int i = pIBegin; i <= pIEnd; ++i) {
            for (int j = pJBegin; j <= pJEnd; ++j) {
                discretization_->p(i,j) += alpha_ * direction_(i,j); // p += alpha * p_k (direction_)
                residual_(i,j)         -= alpha_ * w_(i,j);           // r -= alpha * w
            }
        }

        // apply preconditioner to new residual: compute z_{k+1} into direction_
        applyBlockSGSPreconditioner(num_sgs_sweeps);

        // compute residualNew2_ = r^T z (local) and reduce (wait before other comms)
        residualNew2_ = 0.0;
        for (int i = pIBegin; i <= pIEnd; ++i) {
            for (int j = pJBegin; j <= pJEnd; ++j) {
                residualNew2_ += residual_(i,j) * direction_(i,j); // r^T z_new (local)
            }
        }

        MPI_Request request_residualNew2;
        MPI_Iallreduce(MPI_IN_PLACE, &residualNew2_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_residualNew2);
        MPI_Wait(&request_residualNew2, MPI_STATUS_IGNORE);
        // now safe to communicate halos for direction
        communicateAndSetBoundaryValuesForDirection();

        if (!std::isfinite(residualNew2_)) {
            if (world_rank == 0) std::cerr << "ParallelCG::solve - nonfinite residualNew2_: " << residualNew2_ << "\n";
            break;
        }

        // compute beta and update search direction: p_{k+1} = z_{k+1} + beta * p_k
        double beta = 0.0;
        if (residualOld2_ != 0.0 && std::isfinite(residualOld2_)) {
            beta = residualNew2_ / residualOld2_;
        } else {
            beta = 0.0;
        }
        if (!std::isfinite(beta)) {
            if (world_rank == 0) std::cerr << "ParallelCG::solve - nonfinite beta: " << beta << "\n";
            break;
        }

        // update direction_ in-place using dir_old (which stored p_k)
        idx = 0;
        for (int i = pIBegin; i <= pIEnd; ++i) {
            for (int j = pJBegin; j <= pJEnd; ++j) {
                direction_(i,j) = direction_(i,j) + beta * dir_old[idx++];
            }
        }
        // sanity
        assert(idx == localSize && "dir_old consumption mismatch");

        // finalize for next iteration
        residualOld2_ = residualNew2_;

        // guard against NaN/Inf in residualOld2_
        if (!std::isfinite(residualOld2_)) {
            if (world_rank == 0) std::cerr << "ParallelCG::solve - residualOld2_ became nonfinite at iter " << iter << "\n";
            break;
        }
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