#include "domainCG.hpp"
#include "domain/domain.hpp"
#include <mpi.h>

DomainCG::DomainCG(std::shared_ptr<Discretization> discretization, double epsilon,int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning, std::shared_ptr<Domain> domain) :
    DomainPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning, domain),
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
    {
        const int pIBegin = discretization_->pIBegin();
        const int pIEnd = discretization_->pIEnd();
        const int pJBegin = discretization_->pJBegin();
        const int pJEnd = discretization_->pJEnd();
        sendBufferTopDirection_ = std::vector<double>(pIEnd - pIBegin + 1, 0.0);
        sendBufferBottomDirection_ = std::vector<double>(pIEnd - pIBegin + 1, 0.0);
        sendBufferLeftDirection_ = std::vector<double>(pJEnd - pJBegin + 1, 0.0);
        sendBufferRightDirection_ = std::vector<double>(pJEnd - pJBegin + 1, 0.0);
        recvBufferTopDirection_ = std::vector<double>(pIEnd - pIBegin + 1, 0.0);
        recvBufferBottomDirection_ = std::vector<double>(pIEnd - pIBegin + 1, 0.0);
        recvBufferLeftDirection_ = std::vector<double>(pJEnd - pJBegin + 1, 0.0);
        recvBufferRightDirection_ = std::vector<double>(pJEnd - pJBegin + 1, 0.0);
    }

void DomainCG::solve() {
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
    std::vector<CellInfo> fluidCellsInfo = domain_->getInfoListFluid();

    // check compatibility condition for Neumann Poisson problem and enforce zero-mean RHS
    double bsum_local = 0.0;
    long nfluid_local = 0;

    for (const CellInfo& cellInfo : fluidCellsInfo) {
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];
        bsum_local += discretization_->rhs(i,j);
        nfluid_local++;
    }

    double bsum_global = 0.0;
    long nfluid_global = 0;

    MPI_Allreduce(&bsum_local, &bsum_global, 1, MPI_DOUBLE, MPI_SUM, cartComm_);
    MPI_Allreduce(&nfluid_local, &nfluid_global, 1, MPI_LONG,   MPI_SUM, cartComm_);

    const double bmean = bsum_global / static_cast<double>(nfluid_global);

    // shift RHS so that its mean over all fluid cells is zero: sum(b) = 0
    for (const CellInfo& cellInfo : fluidCellsInfo) {
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];
        discretization_->rhs(i,j) -= bmean;
    }

    // ensure pressure ghost values are up to date before computing A*p
    communicateGhostValues();

    for (const CellInfo& cellInfo : fluidCellsInfo) {
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];

        const double p_i_j = discretization_->p(i,j);
        double p_ip1_j = discretization_->p(i+1,j);
        double p_im1_j = discretization_->p(i-1,j);
        double p_i_jp1 = discretization_->p(i,j+1);
        double p_i_jm1 = discretization_->p(i,j-1);
        // handle boundary faces (Neumann BCs)
        if (cellInfo.faceTop.isBoundaryFace()) {
            p_i_jp1 = p_i_j;
        }
        if (cellInfo.faceBottom.isBoundaryFace()) {
            p_i_jm1 = p_i_j;
        }
        if (cellInfo.faceLeft.isBoundaryFace()) {
            p_im1_j = p_i_j;
        }
        if (cellInfo.faceRight.isBoundaryFace()) {
            p_ip1_j = p_i_j;
        }

        residual_(i,j) = discretization_->rhs(i,j)
            - ((p_ip1_j - 2.0 * p_i_j + p_im1_j) / dx2_
             + (p_i_jp1 - 2.0 * p_i_j + p_i_jm1) / dy2_);
        
        direction_(i,j) = residual_(i,j) * diag_precond; // apply scalar preconditioner: z = M^{-1} r
        residualOld2_ += residual_(i,j) * direction_(i,j); // r^T z

    }

    // global sum of initial r^T z
    MPI_Request request_residual;
    MPI_Iallreduce(MPI_IN_PLACE, &residualOld2_, 1, MPI_DOUBLE, MPI_SUM, cartComm_, &request_residual);
    // need direction halos for first A*direction product later
    communicateAndSetBoundaryValuesForDirection(); //TODOO
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
        for (const CellInfo& cellInfo : fluidCellsInfo) {
            int i = cellInfo.cellIndexPartition[0];
            int j = cellInfo.cellIndexPartition[1];

            const double dir_ij = direction_(i,j);
            double dir_ip1_j = direction_(i+1,j);
            double dir_im1_j = direction_(i-1,j);
            double dir_i_jp1 = direction_(i,j+1);
            double dir_i_jm1 = direction_(i,j-1);

            // handle boundary faces (Neumann BCs)
            if (cellInfo.faceTop.isBoundaryFace()) {
                dir_i_jp1 = dir_ij;
            }
            if (cellInfo.faceBottom.isBoundaryFace()) {
                dir_i_jm1 = dir_ij;
            }
            if (cellInfo.faceLeft.isBoundaryFace()) {
                dir_im1_j = dir_ij;
            }
            if (cellInfo.faceRight.isBoundaryFace()) {
                dir_ip1_j = dir_ij;
            }
            const double w_ij = ((dir_ip1_j - 2.0 * dir_ij + dir_im1_j) / dx2_
                              + (dir_i_jp1 - 2.0 * dir_ij + dir_i_jm1) / dy2_);
            w_(i,j) = w_ij;
            local_dots[0] += dir_ij * w_ij; // dTw_local = z^T w
            local_dots[1] += w_ij * w_ij;   // w
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
        for (const CellInfo& cellInfo : fluidCellsInfo) {
            int i = cellInfo.cellIndexPartition[0];
            int j = cellInfo.cellIndexPartition[1];

            discretization_->p(i,j) += alpha_ * direction_(i,j);
            residual_(i,j) -= alpha_ * w_(i,j);
            const double temp_ij = residual_(i,j) * diag_precond; // z_new = M^{-1} r_new
            direction_(i,j) = temp_ij + beta * direction_(i,j);
        }
        

        // exchange halos for updated direction before next A*direction
        communicateAndSetBoundaryValuesForDirection();

        // update residualOld2_ for next iteration
        residualOld2_ = residualNew2_;
    }

    // final halo exchange for p (solution)
    communicateGhostValues();
    zeroMeanPressure();
    
}



void DomainCG::communicateAndSetBoundaryValuesForDirection() {
    //TODO: Idea for improvement: optimize setting of boundary values by only sending/receiving the necessary values instead of the whole rows/columns (only red boundary values needed or black)
    
    const int pIBegin = discretization_->pIBegin();
    const int pIEnd = discretization_->pIEnd();
    const int pJBegin = discretization_->pJBegin();
    const int pJEnd = discretization_->pJEnd();

    // init MPI request variables
    MPI_Request requestTop, requestBottom, requestLeft, requestRight;
    MPI_Request recvRequestTop, recvRequestBottom, recvRequestLeft, recvRequestRight;

    // Fill easy boundary conditions that do not need communication
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            direction_(i, pJEnd + 1) = direction_(i, pJEnd);
        }
    } else {
        for (int i = pIBegin; i <= pIEnd; i++) {
            sendBufferTopDirection_[i - pIBegin] = direction_(i, pJEnd);
        }
        // instantiate non-blocking sends and receives
        MPI_Isend(sendBufferTopDirection_.data(), sendBufferTopDirection_.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, cartComm_, &requestTop);
        MPI_Irecv(recvBufferTopDirection_.data(), recvBufferTopDirection_.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, cartComm_, &recvRequestTop);
    }

    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            direction_(i, pJBegin - 1) = direction_(i, pJBegin);
        }
    } else {
        for (int i = pIBegin; i <= pIEnd; i++) {
            sendBufferBottomDirection_[i - pIBegin] = direction_(i, pJBegin);
        }
        MPI_Isend(sendBufferBottomDirection_.data(), sendBufferBottomDirection_.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, cartComm_, &requestBottom);
        MPI_Irecv(recvBufferBottomDirection_.data(), recvBufferBottomDirection_.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, cartComm_, &recvRequestBottom);
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = pJBegin - 1; j <= pJEnd + 1; j++) {
            direction_(pIBegin - 1, j) = direction_(pIBegin, j);
        }
    } else {
        for (int j = pJBegin; j <= pJEnd; j++) {
            sendBufferLeftDirection_[j - pJBegin] = direction_(pIBegin, j);
        }
        MPI_Isend(sendBufferLeftDirection_.data(), sendBufferLeftDirection_.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, cartComm_, &requestLeft);
        MPI_Irecv(recvBufferLeftDirection_.data(), recvBufferLeftDirection_.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, cartComm_, &recvRequestLeft);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = pJBegin - 1; j <= pJEnd + 1; j++) {
            direction_(pIEnd + 1, j) = direction_(pIEnd, j);
        }
    } else {
        for (int j = pJBegin; j <= pJEnd; j++) {
            sendBufferRightDirection_[j - pJBegin] = direction_(pIEnd, j);
        }
        MPI_Isend(sendBufferRightDirection_.data(), sendBufferRightDirection_.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, cartComm_, &requestRight);
        MPI_Irecv(recvBufferRightDirection_.data(), recvBufferRightDirection_.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, cartComm_, &recvRequestRight);
    }

    // nachdem kommuniziert wurde setzren wir die ghost nodes mit den empfangenen werten für ränder die nicht am globalen rand liegen
    
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&recvRequestTop, MPI_STATUS_IGNORE);  // Wait for receive
        for (int i = pIBegin; i <= pIEnd; i++) {
            // set ghost cells
            direction_(i, pJEnd + 1) = recvBufferTopDirection_[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&recvRequestBottom, MPI_STATUS_IGNORE);  // Wait for receive
        for (int i = pIBegin; i <= pIEnd; i++) {
            // set ghost cells
            direction_(i, pJBegin - 1) = recvBufferBottomDirection_[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&recvRequestLeft, MPI_STATUS_IGNORE);  // Wait for receive
        for (int j = pJBegin; j <= pJEnd; j++) {
            // set ghost cells
            direction_(pIBegin - 1, j) = recvBufferLeftDirection_[j - pJBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&recvRequestRight, MPI_STATUS_IGNORE);  // Wait for receive
        for (int j = pJBegin; j <= pJEnd; j++) {
            // set ghost cells
            direction_(pIEnd + 1, j) = recvBufferRightDirection_[j - pJBegin];
        }
    }
}