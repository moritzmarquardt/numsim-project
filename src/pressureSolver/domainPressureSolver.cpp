#include "pressureSolver/domainPressureSolver.hpp"
#include "pressureSolver/PressureSolver.hpp"
#include "discretization/discretization.hpp"
#include "partitioning/partitioning.hpp"
#include <cmath>

DomainPressureSolver::DomainPressureSolver(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning, std::shared_ptr<Domain> domain) :
    PressureSolver(discretization, epsilon, maximumNumberOfIterations), partitioning_(partitioning), domain_(domain) {
    cartComm_ = partitioning_->getCartComm();
}

void DomainPressureSolver::computeResidualNorm() {
    double residual_norm_squared = 0.0;
    const int N_global = partitioning_->getNCellsGlobal(); // global number of real cells
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();

    std::vector<CellInfo> fluidCellsInfo = domain_->getInfoListFluid();
    int n = fluidCellsInfo.size();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = fluidCellsInfo[idx];
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];

        double p_i_j = discretization_->p(i,j);
        double p_ip1_j = discretization_->p(i+1,j);
        double p_im1_j = discretization_->p(i-1,j);
        double p_i_jp1 = discretization_->p(i,j+1);
        double p_i_jm1 = discretization_->p(i,j-1);

        if (cellInfo.topIsBoundaryFace()) {
            p_i_jp1 = p_i_j;
        }
        if (cellInfo.bottomIsBoundaryFace()) {
            p_i_jm1 = p_i_j;
        }
        if (cellInfo.leftIsBoundaryFace()) {
            p_im1_j = p_i_j;
        }
        if (cellInfo.rightIsBoundaryFace()) {
            p_ip1_j = p_i_j;
        }

        const double rhs_ij = discretization_->rhs(i,j);
        const double p_xx = (p_ip1_j - 2.0 * p_i_j + p_im1_j) / dx2;
        const double p_yy = (p_i_jp1 - 2.0 * p_i_j + p_i_jm1) / dy2;
        const double residual_ij = rhs_ij - (p_xx + p_yy);
        residual_norm_squared += residual_ij * residual_ij;
    }

    double residual_norm_global = 0.0;
    // perform global reduction to sum up residual norms from all processes
    // blocking call is acceptable here since we need the result before proceeding
    MPI_Allreduce(&residual_norm_squared, &residual_norm_global, 1, MPI_DOUBLE, MPI_SUM, cartComm_);
    residualNorm_ = residual_norm_global / N_global;
}

void DomainPressureSolver::communicateGhostValues() {

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

    std::vector<double> recvBufferTop(pIEnd - pIBegin + 1, 0.0);
    std::vector<double> recvBufferBottom(pIEnd - pIBegin + 1, 0.0);
    std::vector<double> recvBufferLeft(pJEnd - pJBegin + 1, 0.0);
    std::vector<double> recvBufferRight(pJEnd - pJBegin + 1, 0.0);

    // init MPI request variables
    MPI_Request requestTop, requestBottom, requestLeft, requestRight;
    MPI_Request recvRequestTop, recvRequestBottom, recvRequestLeft, recvRequestRight;

    // Fill easy boundary conditions that do not need communication
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            sendBufferTop[i - pIBegin] = discretization_->p(i, pJEnd);
        }
        // instantiate non-blocking sends and receives
        MPI_Isend(sendBufferTop.data(), sendBufferTop.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, cartComm_, &requestTop);
        MPI_Irecv(recvBufferTop.data(), recvBufferTop.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, cartComm_, &recvRequestTop); 
        // override request top with the receive request. and it does not matter since we only wait for receives later (if a receive is not done yet, the send cannot be done either)
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            sendBufferBottom[i - pIBegin] = discretization_->p(i, pJBegin);
        }
        MPI_Isend(sendBufferBottom.data(), sendBufferBottom.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, cartComm_, &requestBottom);
        MPI_Irecv(recvBufferBottom.data(), recvBufferBottom.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, cartComm_, &recvRequestBottom);
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = pJBegin; j <= pJEnd; j++) {
            sendBufferLeft[j - pJBegin] = discretization_->p(pIBegin, j);
        }
        MPI_Isend(sendBufferLeft.data(), sendBufferLeft.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, cartComm_, &requestLeft);
        MPI_Irecv(recvBufferLeft.data(), recvBufferLeft.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, cartComm_, &recvRequestLeft);
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = pJBegin; j <= pJEnd; j++) {
            sendBufferRight[j - pJBegin] = discretization_->p(pIEnd, j);
        }
        MPI_Isend(sendBufferRight.data(), sendBufferRight.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, cartComm_, &requestRight);
        MPI_Irecv(recvBufferRight.data(), recvBufferRight.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, cartComm_, &recvRequestRight);
    }

    // nachdem kommuniziert wurde setzren wir die ghost nodes mit den empfangenen werten für ränder die nicht am globalen rand liegen
    
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&recvRequestTop, MPI_STATUS_IGNORE);  // Wait for receive
        for (int i = pIBegin; i <= pIEnd; i++) {
            // set ghost cells
            discretization_->p(i, pJEnd + 1) = recvBufferTop[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&recvRequestBottom, MPI_STATUS_IGNORE);  // Wait for receive
        for (int i = pIBegin; i <= pIEnd; i++) {
            // set ghost cells
            discretization_->p(i, pJBegin - 1) = recvBufferBottom[i - pIBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&recvRequestLeft, MPI_STATUS_IGNORE);  // Wait for receive
        for (int j = pJBegin; j <= pJEnd; j++) {
            // set ghost cells
            discretization_->p(pIBegin - 1, j) = recvBufferLeft[j - pJBegin];
        }
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&recvRequestRight, MPI_STATUS_IGNORE);  // Wait for receive
        for (int j = pJBegin; j <= pJEnd; j++) {
            // set ghost cells
            discretization_->p(pIEnd + 1, j) = recvBufferRight[j - pJBegin];
        }
    }
}