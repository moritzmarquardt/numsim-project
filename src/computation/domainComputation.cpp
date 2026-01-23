#include "domainComputation.hpp"
#include <cmath>

void DomainComputation::initialize(int argc, char *argv[]) {
    
    settings_.loadFromFile(argv[1]);

    partitioning_ = std::make_shared<Partitioning>();
    partitioning_->initialize(settings_.nCells);
    std::array<int,2> nCellsLocal = partitioning_->nCellsLocal();

    // calculate mesh width
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    domain_ = Domain(&settings_, partitioning_);
    domain_.readDomainFile(argv[2]);
    // create discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(nCellsLocal, meshWidth_, settings_.alpha, partitioning_);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(nCellsLocal, meshWidth_,partitioning_);
    }

    // create pressure solver
    // TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<RedBlackGaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_, domain_);
    } else if (settings_.pressureSolver == "SOR") {
        // TODO: calculate optimal omega (not in settings hardcoded)
        pressureSolver_ = std::make_unique<RedBlackSOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega, partitioning_);
    } else if (settings_.pressureSolver == "CG") {
        pressureSolver_ = std::make_unique<ParallelCG>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);
    } else {
        std::cerr << "Error: Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(EXIT_FAILURE);
    }
    // // Hard code CG
    // pressureSolver_ = std::make_unique<ParallelCG>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);

    // create output writers
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, *partitioning_);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, *partitioning_);

    //Set cartComm_ member
    cartComm_ = partitioning_->getCartComm();
}

void DomainComputation::runSimulation() {
    applyInitialBoundaryValues();

    double currentTime = 0.0;
    int iterationCount = 0;
    const double time_eps = 1e-8;
    int nOutputs = 1;

    while (currentTime < settings_.endTime - time_eps) {
        communicateGhostCells();

        computeTimeStepWidth();

        if (currentTime + dt_ > settings_.endTime - time_eps) {
            dt_ = settings_.endTime - currentTime;
        }


        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();

        currentTime += dt_;
        iterationCount++;

        // this was the fix!!!
        if (currentTime >= nOutputs) {
            outputWriterParaview_->writeFile(currentTime);
            nOutputs = nOutputs + 1;
        }       
    }
}

void DomainComputation::computeTimeStepWidth() {
    // compute time step width based on CFL condition
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    double dx2 = dx * dx;
    double dy2 = dy * dy;  

    double dt_diff_cond = 0.5 * settings_.re * (dx2 * dy2) / (dx2 + dy2);
    double dt_conv_cond_local = std::min(dx / discretization_->u().computeMaxAbs(), dy / discretization_->v().computeMaxAbs());

    MPI_Request time_request;
    double dt_conv_cond_global = 0.0;
    // perform global reduction to find minimum dt_conv_cond across all processes
    MPI_Iallreduce(&dt_conv_cond_local, &dt_conv_cond_global, 1, MPI_DOUBLE, MPI_MIN, cartComm_, &time_request);
    MPI_Wait(&time_request, MPI_STATUS_IGNORE);

    double dt_prelim = settings_.tau * std::min(dt_diff_cond, dt_conv_cond_global);

    if (dt_prelim < settings_.maximumDt) {
        dt_ = dt_prelim;
    } else {
        dt_ = settings_.maximumDt;
        std::cout << "Warning: Time step width limited by maximumDt!" << std::endl;
    }
}

void DomainComputation::applyInitialBoundaryValues() {
    //it is sufficient to only go though all cells and only look at the right and top face since then we will go through all faces exactly once. the values set here are only when we have dirichlet BCs directly orthogonally flowing through the face direction. These are onyl set once in the beginning and then never touched again.
    std::vector<CellInfo> allCellsInfo = domain_.getInfoListAll();
    int n = allCellsInfo.size();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = allCellsInfo[idx];
        if (cellInfo.rightIsBoundaryFace() || cellInfo.topIsBoundaryFace()) {
            int i = cellInfo.cellIndexPartition[0];
            int j = cellInfo.cellIndexPartition[1];
            double faceBCValue; // only stores the value of the faceBC that is orthogonal to face orientation.
            // top face
            if (cellInfo.topIsBoundaryFace()) {
                faceBCValue = cellInfo.faceTop.dirichletV.value();
                if (faceBCValue != 0.0) {
                    discretization_->v(i, j) = faceBCValue;
                    discretization_->g(i, j) = faceBCValue;
                }
            }
            // right face
            if (cellInfo.rightIsBoundaryFace()) {
                faceBCValue = cellInfo.faceRight.dirichletU.value();
                if (faceBCValue != 0.0) {
                    discretization_->u(i, j) = faceBCValue;
                    discretization_->f(i, j) = faceBCValue;
                }
            }
        }
    }
}


void DomainComputation::communicateGhostCells() {
    const int uIBegin = discretization_->uIBegin();
    const int uIEnd = discretization_->uIEnd();
    const int uJBegin = discretization_->uJBegin();
    const int uJEnd = discretization_->uJEnd();

    const int vIBegin = discretization_->vIBegin();
    const int vIEnd = discretization_->vIEnd();
    const int vJBegin = discretization_->vJBegin();
    const int vJEnd = discretization_->vJEnd();

    // buffers for sending and receiving data
    std::vector<double> sendBufferTopU(uIEnd - uIBegin + 1, 0.0);
    std::vector<double> sendBufferTopV(vIEnd - vIBegin + 1, 0.0);
    std::vector<double> sendBufferBottomU(uIEnd - uIBegin + 1, 0.0);
    std::vector<double> sendBufferBottomV(vIEnd - vIBegin + 1, 0.0);
    std::vector<double> sendBufferLeftU(uJEnd - uJBegin + 1, 0.0);
    std::vector<double> sendBufferLeftV(vJEnd - vJBegin + 1, 0.0);
    std::vector<double> sendBufferRightU(uJEnd - uJBegin + 1, 0.0);
    std::vector<double> sendBufferRightV(vJEnd - vJBegin + 1, 0.0);

    MPI_Request requestsTop, requestsBottom, requestsLeft, requestsRight;
    const int TAG_U = 0;
    const int TAG_V = 1;

    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        // otherwise communicate with the top neighbour
        for (int i = uIBegin; i <= uIEnd; i++) {
            sendBufferTopU[i - uIBegin] = discretization_->u(i,uJEnd);
        }

        for (int i = vIBegin; i <= vIEnd; i++) {
            sendBufferTopV[i - vIBegin] = discretization_->v(i,vJEnd - 1);
        }
        // instantiate non-blocking sends and receives
        MPI_Isend(sendBufferTopU.data(), sendBufferTopU.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), TAG_U, cartComm_, &requestsTop);
        MPI_Isend(sendBufferTopV.data(), sendBufferTopV.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), TAG_V, cartComm_, &requestsTop);

        MPI_Irecv(sendBufferTopU.data(), sendBufferTopU.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), TAG_U, cartComm_, &requestsTop);
        MPI_Irecv(sendBufferTopV.data(), sendBufferTopV.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), TAG_V, cartComm_, &requestsTop);
    }
    
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = uIBegin; i <= uIEnd; i++) {
            sendBufferBottomU[i - uIBegin] = discretization_->u(i,uJBegin);
        }

        for (int i = vIBegin; i <= vIEnd; i++) {
            sendBufferBottomV[i - vIBegin] = discretization_->v(i,vJBegin + 1); // +1 because we have two layers of gjost cells at the bottom (just like at the left)
        }
        MPI_Isend(sendBufferBottomU.data(), sendBufferBottomU.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_U, cartComm_, &requestsBottom);
        MPI_Isend(sendBufferBottomV.data(), sendBufferBottomV.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_V, cartComm_, &requestsBottom);

        MPI_Irecv(sendBufferBottomU.data(), sendBufferBottomU.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_U, cartComm_, &requestsBottom);
        MPI_Irecv(sendBufferBottomV.data(), sendBufferBottomV.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_V, cartComm_, &requestsBottom);
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = uJBegin; j <= uJEnd; j++) {
            sendBufferLeftU[j - uJBegin] = discretization_->u(uIBegin + 1,j); // +1 because we have two layers of ghost cells at the left (just like at the bottom)
        }

        for (int j = vJBegin; j <= vJEnd; j++) {
            sendBufferLeftV[j - vJBegin] = discretization_->v(vIBegin,j);
        }
        MPI_Isend(sendBufferLeftU.data(), sendBufferLeftU.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), TAG_U, cartComm_, &requestsLeft);
        MPI_Isend(sendBufferLeftV.data(), sendBufferLeftV.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), TAG_V, cartComm_, &requestsLeft);
        
        MPI_Irecv(sendBufferLeftU.data(), sendBufferLeftU.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), TAG_U, cartComm_, &requestsLeft);
        MPI_Irecv(sendBufferLeftV.data(), sendBufferLeftV.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), TAG_V, cartComm_, &requestsLeft);
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = uJBegin; j <= uJEnd; j++) {
            sendBufferRightU[j - uJBegin] = discretization_->u(uIEnd - 1,j);
        }

        for (int j = vJBegin; j <= vJEnd; j++) {
            sendBufferRightV[j - vJBegin] = discretization_->v(vIEnd,j);
        }
        MPI_Isend(sendBufferRightU.data(), sendBufferRightU.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), TAG_U, cartComm_, &requestsRight);
        MPI_Isend(sendBufferRightV.data(), sendBufferRightV.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), TAG_V, cartComm_, &requestsRight);

        MPI_Irecv(sendBufferRightU.data(), sendBufferRightU.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), TAG_U, cartComm_, &requestsRight);
        MPI_Irecv(sendBufferRightV.data(), sendBufferRightV.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), TAG_V, cartComm_, &requestsRight);
    }

    // wait for all communications to finish and set ghost values
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestsTop, MPI_STATUS_IGNORE);
        for (int i = uIBegin; i <= uIEnd; i++) {
            discretization_->u(i,uJEnd + 1) = sendBufferTopU[i - uIBegin];
        }
        for (int i = vIBegin; i <= vIEnd; i++) {
            discretization_->v(i,vJEnd + 1) = sendBufferTopV[i - vIBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestsBottom, MPI_STATUS_IGNORE);
        for (int i = uIBegin; i <= uIEnd; i++) {
            discretization_->u(i,uJBegin - 1) = sendBufferBottomU[i - uIBegin];
        }
        for (int i = vIBegin; i <= vIEnd; i++) {
            discretization_->v(i,vJBegin - 1) = sendBufferBottomV[i - vIBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestsLeft, MPI_STATUS_IGNORE);
        for (int j = uJBegin; j <= uJEnd; j++) {
            discretization_->u(uIBegin - 1,j) = sendBufferLeftU[j - uJBegin];
        }
        for (int j = vJBegin; j <= vJEnd; j++) {
            discretization_->v(vIBegin - 1,j) = sendBufferLeftV[j - vJBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestsRight, MPI_STATUS_IGNORE);
        for (int j = uJBegin; j <= uJEnd; j++) {
            discretization_->u(uIEnd + 1,j) = sendBufferRightU[j - uJBegin];
        }
        for (int j = vJBegin; j <= vJEnd; j++) {
            discretization_->v(vIEnd + 1,j) = sendBufferRightV[j - vJBegin];
        }
    }
}


double DomainComputation::computeD2uDx2(double u_ip1_j, double u_i_j, double u_im1_j) const {
    const double dx = discretization_->dx();
    return (u_ip1_j - 2.0 * u_i_j + u_im1_j) / (dx * dx);
}

double DomainComputation::computeD2uDy2(double u_i_jp1, double u_i_j, double u_i_jm1) const {
    const double dy = discretization_->dy();
    return (u_i_jp1 - 2.0 * u_i_j + u_i_jm1) / (dy * dy);
}

double DomainComputation::computeD2vDx2(double v_ip1_j, double v_i_j, double v_im1_j) const {
    const double dx = discretization_->dx();
    return (v_ip1_j - 2.0 * v_i_j + v_im1_j) / (dx * dx);
}

double DomainComputation::computeD2vDy2(double v_i_jp1, double v_i_j, double v_i_jm1) const {
    const double dy = discretization_->dy();
    return (v_i_jp1 - 2.0 * v_i_j + v_i_jm1) / (dy * dy);
}

double DomainComputation::computeDpDx(double p_ip1_j, double p_i_j) const {
    const double dx = discretization_->dx();
    return (p_ip1_j - p_i_j) / dx;
}

double DomainComputation::computeDpDy(double p_i_jp1, double p_i_j) const {
    const double dy = discretization_->dy();
    return (p_i_jp1 - p_i_j) / dy;
}

double DomainComputation::computeDu2Dx(double u_i_j, double u_im1_j, double u_ip1_j) const {
    const double u_right_sum = (u_i_j + u_ip1_j) / 2.0;
    const double u_left_sum = (u_im1_j + u_i_j) / 2.0;
    const double u_right_diff = (u_i_j - u_ip1_j) / 2.0;
    const double u_left_diff = (u_im1_j - u_i_j) / 2.0;
    const double dx = discretization_->dx();

    const double central_diff_term = ((u_right_sum*u_right_sum) - (u_left_sum*u_left_sum)) / dx;
    const double donor_cell_term = (std::fabs(u_right_sum) * u_right_diff - std::fabs(u_left_sum) * u_left_diff) / dx;
    
    return central_diff_term + settings_.alpha * donor_cell_term;
}

double DomainComputation::computeDv2Dy(double v_i_j, double v_i_jp1, double v_i_jm1) const {
    const double v_top_sum = (v_i_j + v_i_jp1) / 2.0;
    const double v_bottom_sum = (v_i_jm1 + v_i_j) / 2.0;
    const double v_top_diff = (v_i_j - v_i_jp1) / 2.0;
    const double v_bottom_diff = (v_i_jm1 - v_i_j) / 2.0;
    const double dy = discretization_->dy();

    const double central_diff_term = ((v_top_sum*v_top_sum) - (v_bottom_sum*v_bottom_sum)) / dy;
    const double donor_cell_term = (std::fabs(v_top_sum) * v_top_diff - std::fabs(v_bottom_sum) * v_bottom_diff) / dy;

    return central_diff_term + settings_.alpha * donor_cell_term;
}

double DomainComputation::computeDuvDx(double u_i_j, double u_i_jp1, double u_im1_j, double u_im1_jp1, double v_i_j, double v_ip1_j, double v_im1_j) const {
    const double u_top_sum = (u_i_j + u_i_jp1) / 2.0;
    const double v_right_sum = (v_i_j + v_ip1_j) / 2.0;
    const double u_top_left_sum = (u_im1_j + u_im1_jp1) / 2.0;
    const double v_left_sum = (v_im1_j + v_i_j) / 2.0;
    const double v_right_diff = (v_i_j - v_ip1_j) / 2.0;
    const double v_left_diff = (v_im1_j - v_i_j) / 2.0;
    const double dx = discretization_->dx();

    const double central_diff_term = (u_top_sum*v_right_sum - u_top_left_sum*v_left_sum) / dx;
    const double donor_cell_term = (std::fabs(u_top_sum) * v_right_diff - std::fabs(u_top_left_sum) * v_left_diff) / dx;

    return central_diff_term + settings_.alpha * donor_cell_term;
}

double DomainComputation::computeDuvDy(double u_i_j, double u_i_jp1, double u_i_jm1, double v_i_j, double v_ip1_j, double v_i_jm1, double v_ip1_jm1) const {
    const double v_top_right_sum = (v_i_j + v_ip1_j) / 2.0;
    const double u_top_right_sum = (u_i_j + u_i_jp1) / 2.0;
    const double v_right_bottom_sum = (v_i_jm1 + v_ip1_jm1) / 2.0;
    const double u_bottom_right_sum = (u_i_jm1 + u_i_j) / 2.0;
    const double u_top_diff = (u_i_j - u_i_jp1) / 2.0;
    const double u_bottom_diff = (u_i_jm1 - u_i_j) / 2.0;
    const double dy = discretization_->dy();

    const double central_diff_term = (v_top_right_sum*u_top_right_sum - v_right_bottom_sum*u_bottom_right_sum) / dy;
    const double donor_cell_term = (std::fabs(v_top_right_sum) * u_top_diff - std::fabs(v_right_bottom_sum) * u_bottom_diff) / dy;

    return central_diff_term + settings_.alpha * donor_cell_term;
}


void DomainComputation::computePreliminaryVelocities() {
    std::vector<CellInfo> allCellsInfo = domain_.getInfoListFluid();
    int n = allCellsInfo.size();
    double dx = discretization_->dx();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = allCellsInfo[idx];
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];
        if (cellInfo.hasAnyBoundaryFace()) {
            // check which faces are boundary faces and apply one sided stencils accordingly
            bool topBC = cellInfo.topIsBoundaryFace();
            bool rightBC = cellInfo.rightIsBoundaryFace();
            bool bottomBC = cellInfo.bottomIsBoundaryFace();
            bool leftBC = cellInfo.leftIsBoundaryFaceC();

            double u_i_j = discretization_->u(i,j);
            double u_ip1_j = discretization_->u(i+1,j);
            double u_im1_j = discretization_->u(i-1,j);
            double u_i_jp1 = discretization_->u(i,j+1);
            double u_i_jm1 = discretization_->u(i,j-1);
            double u_im1_jp1 = discretization_->u(i-1,j+1);
            double v_i_j = discretization_->v(i,j);
            double v_ip1_j = discretization_->v(i+1,j);
            double v_im1_j = discretization_->v(i-1,j);
            double v_i_jm1 = discretization_->v(i,j-1);
            double v_ip1_jm1 = discretization_->v(i+1,j-1);
            double v_i_jp1 = discretization_->v(i,j+1);
            bool calcB = true;
            bool calcA = true;
            double dy = discretization_->dy();
            double dx = discretization_->dx();

            if (topBC) {
                if (cellInfo.faceTop.dirichletU.has_value()) {
                    u_i_jp1 = 2.0 * cellInfo.faceTop.dirichletU.value() - u_i_j;
                }  else if (cellInfo.faceTop.neumannU.has_value()) {
                    u_i_jp1 = u_i_j + cellInfo.faceTop.neumannU.value() * dy;                    
                }
                if (cellInfo.faceTop.dirichletV.has_value()) {
                    calcB = false;
                } else if (cellInfo.faceTop.neumannV.has_value()) {
                    discretization_->g(i,j) = v_i_jm1 + cellInfo.faceTop.neumannV.value() * dy;
                    calcB = false;
                }
            }
            if (rightBC) {
                if (cellInfo.faceRight.dirichletU.has_value()) {
                    calcA = false;
                } else if (cellInfo.faceRight.neumannU.has_value()) {
                    discretization_->f(i,j) = u_im1_j + cellInfo.faceRight.neumannU.value() * dx;
                    calcA = false;                 

                }
                if (cellInfo.faceRight.dirichletV.has_value()) {
                    v_ip1_j = 2 * cellInfo.faceRight.dirichletV.value() - v_i_j;
                }  else if (cellInfo.faceRight.neumannV.has_value()) {
                    v_ip1_j = v_i_j + cellInfo.faceRight.neumannV.value() * dx;
                }
            }
            if (leftBC) {
                if (cellInfo.faceLeft.dirichletU.has_value()) {
                    // do nothing
                    continue;
                } else if (cellInfo.faceLeft.neumannU.has_value()) { //neuman gewixe
                    u_im1_j = u_i_j + cellInfo.faceLeft.neumannU.value() * dx;
                    discretization_->f(i-1,j) = u_im1_j; // set f to the neumann value (corresponds to a solid cell bc neumann left means solid obstacle to the left)
                } 
                if (cellInfo.faceLeft.dirichletV.has_value()) {
                    v_im1_j = 2.0 * cellInfo.faceLeft.dirichletV.value() - v_i_j;
                }  else if (cellInfo.faceLeft.neumannV.has_value()) {
                    v_im1_j = v_i_j + cellInfo.faceLeft.neumannV.value() * dx;
                }
            }
            if (bottomBC) {
                if (cellInfo.faceBottom.dirichletU.has_value()) {
                    u_i_jm1 = 2.0 * cellInfo.faceBottom.dirichletU.value() - u_i_j;
                }  else if (cellInfo.faceBottom.neumannU.has_value()) {
                    u_i_jm1 = u_i_j + cellInfo.faceBottom.neumannU.value() * dy;     
                }
                if (cellInfo.faceBottom.dirichletV.has_value()) {
                    // do nothing
                    continue;
                } else if (cellInfo.faceBottom.neumannV.has_value()) {
                    v_i_jm1 = v_i_j + cellInfo.faceBottom.neumannV.value() * dy;
                    discretization_->g(i,j-1) = v_i_jm1; // set g to the neumann value
                }
            }

            
            if (calcA) {
                double A_ij = 1 / settings_.re * (computeD2uDx2(u_ip1_j, u_i_j, u_im1_j) + computeD2uDy2(u_i_jp1, u_i_j, u_i_jm1)) - computeDu2Dx(u_i_j, u_im1_j, u_ip1_j) - computeDuvDy(u_i_j, u_i_jp1, u_i_jm1, v_i_j, v_ip1_j, v_i_jm1, v_ip1_jm1) + settings_.g[0];
                discretization_->f(i,j) = u_i_j + A_ij * dt_;
            }

            if (calcB) {
                double B_ij = 1 / settings_.re * (computeD2vDx2(v_ip1_j, v_i_j, v_im1_j) + computeD2vDy2(v_i_jp1, v_i_j, v_i_jm1)) - computeDuvDx(u_i_j, u_i_jp1, u_im1_j, u_im1_jp1, v_i_j, v_ip1_j, v_im1_j) - computeDv2Dy(v_i_j, v_i_jp1, v_i_jm1) + settings_.g[1];
                discretization_->g(i,j) = v_i_j + B_ij * dt_;   
            }      

        } else {
            // inner cell, normal stencil
            double A_ij = 1 / settings_.re * ( discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j)) - discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j) + settings_.g[0];
            discretization_->f(i,j) = discretization_->u(i,j) + A_ij * dt_;

            double B_ij = 1 / settings_.re * ( discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j))
            - discretization_->computeDuvDx(i,j) - discretization_->computeDv2Dy(i,j) + settings_.g[1];
            discretization_->g(i,j) = discretization_->v(i,j) + B_ij * dt_;
        }

    }    
}

void DomainComputation::computeRightHandSide() {
    std::vector<CellInfo> allCellsInfo = domain_.getInfoListFluid();
    int n = allCellsInfo.size();
    double dx = discretization_->dx();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = allCellsInfo[idx];
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];

        double rhs_ij = (discretization_->f(i,j) - discretization_->f(i-1,j)) / discretization_->dx() + (discretization_->g(i,j) - discretization_->g(i,j-1)) / discretization_->dy();

        discretization_->rhs(i,j) = rhs_ij / dt_;
    }
}

void DomainComputation::computePressure() {
    pressureSolver_->solve();
}

void DomainComputation::computeVelocities() {
    // update velocities based on new pressure field
    std::vector<CellInfo> fluidCellsInfo = domain_.getInfoListFluid();
    int n = fluidCellsInfo.size();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = fluidCellsInfo[idx];
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];

        double p_i_j = discretization_->p(i,j);
        double p_ip1_j = discretization_->p(i+1,j);
        double p_i_jp1 = discretization_->p(i,j+1);


        if (!cellInfo.rightIsBoundaryFace()){
            discretization_->u(i,j) = discretization_->f(i,j) - dt_ * computeDpDx(p_ip1_j, p_i_j);
        } else if (cellInfo.faceRight.neumannU.has_value()) {
            // right boundary face
            discretization_->u(i,j) = discretization_->f(i,j);
        }

        if (!cellInfo.topIsBoundaryFace()){
            discretization_->v(i,j) = discretization_->g(i,j) - dt_ * computeDpDy(p_i_jp1, p_i_j);
        } else if (cellInfo.faceTop.neumannV.has_value()) {
            // top boundary face
            discretization_->v(i,j) = discretization_->g(i,j);
        }

    }
}

void DomainComputation::computeTimeStepWidth() {
    // compute time step width based on CFL condition
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    double dx2 = dx * dx;
    double dy2 = dy * dy;  

    double dt_diff_cond = 0.5 * settings_.re * (dx2 * dy2) / (dx2 + dy2);
    double dt_conv_cond_local = std::min(dx / discretization_->u().computeMaxAbs(), dy / discretization_->v().computeMaxAbs());

    MPI_Request time_request;
    double dt_conv_cond_global = 0.0;
    // perform global reduction to find minimum dt_conv_cond across all processes
    MPI_Iallreduce(&dt_conv_cond_local, &dt_conv_cond_global, 1, MPI_DOUBLE, MPI_MIN, cartComm_, &time_request);
    MPI_Wait(&time_request, MPI_STATUS_IGNORE);

    double dt_prelim = settings_.tau * std::min(dt_diff_cond, dt_conv_cond_global);

    if (dt_prelim < settings_.maximumDt) {
        dt_ = dt_prelim;
    } else {
        dt_ = settings_.maximumDt;
        std::cout << "Warning: Time step width limited by maximumDt!" << std::endl;
    }
}