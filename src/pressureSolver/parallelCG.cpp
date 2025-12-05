#include "parallelCG.hpp"


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

    rTr_ = 0.0;

    // compute initial residual r = b - A*p and set direction = r
    for (int i = pIBegin; i <= pIEnd; i++) {
        for (int j = pJBegin; j <= pJEnd; j++) {
            residual_(i,j) = discretization_->rhs(i,j) - ((discretization_->p(i+1,j) - 2.0 * discretization_->p(i,j) + discretization_->p(i-1,j)) / dx2_
        + (discretization_->p(i,j+1) - 2.0 * discretization_->p(i,j) + discretization_->p(i,j-1)) / dy2_);
            direction_(i,j) = residual_(i,j);
            rTr_ += residual_(i,j) * direction_(i,j);
        }
    }

    MPI_Request request_residual;
    MPI_Iallreduce(MPI_IN_PLACE, &rTr_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_residual);
    communicateAndSetBoundaryValuesForDirection();
    MPI_Wait(&request_residual, MPI_STATUS_IGNORE);

    if (rTr_ / N < eps_2) {
        return; // initial guess is good enough
    }

    // main CG iteration loop
    double wTw_ = 0.0;
    double rTw_ = 0.0;
    double dTw_ = 0.0;
    while (iter < maximumNumberOfIterations_ && rTr_ / N > eps_2) {
        iter++;

        // compute w = A * direction
        for (int i = pIBegin; i <= pIEnd; i++) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                const double w_ij =  ((direction_(i+1,j) - 2.0 * direction_(i,j) + direction_(i-1,j)) / dx2_
            + (direction_(i,j+1) - 2.0 * direction_(i,j) + direction_(i,j-1)) / dy2_);
                w_(i,j) = w_ij;
                rTw_ += residual_(i,j) * w_ij;
                wTw_ += w_ij * w_ij;
                // rTr_ += residual_(i,j) * residual_(i,j);
                dTw_ += direction_(i,j) * w_ij;
            }
        }

        MPI_Request request_rTw, request_wTw, request_rTr, request_dTw;
        MPI_Iallreduce(MPI_IN_PLACE, &rTw_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_rTw);
        MPI_Iallreduce(MPI_IN_PLACE, &wTw_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_wTw);
        MPI_Iallreduce(MPI_IN_PLACE, &rTr_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_rTr);
        MPI_Iallreduce(MPI_IN_PLACE, &dTw_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_dTw);

        MPI_Wait(&request_rTr, MPI_STATUS_IGNORE);
        MPI_Wait(&request_dTw, MPI_STATUS_IGNORE);
        alpha_ = rTr_ / dTw_;

        for (int i = pIBegin; i <= pIEnd; i++) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                discretization_->p(i,j) += alpha_ * direction_(i,j);
                residual_(i,j) -= alpha_ * w_(i,j);
            }
        }

        MPI_Wait(&request_rTw, MPI_STATUS_IGNORE);
        MPI_Wait(&request_wTw, MPI_STATUS_IGNORE);
        rTrNew_ = rTr_ - 2 * alpha_ * rTw_ + alpha_ * alpha_ * wTw_;

        if (rTrNew_ / N < eps_2) {
            break; // converged
        }

        beta_ = rTrNew_ / rTr_;

        for (int i = pIBegin; i <= pIEnd; i++) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                direction_(i,j) = residual_(i,j) + beta_ * direction_(i,j);
            }
        }
        communicateAndSetBoundaryValuesForDirection();
        rTr_ = rTrNew_;
        wTw_ = 0.0;
        rTw_ = 0.0;
        dTw_ = 0.0;
    }

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
            sendBufferTop[i - pIBegin] = discretization_->p(i, pJEnd);
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