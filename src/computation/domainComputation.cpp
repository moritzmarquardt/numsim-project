#include "domainComputation.hpp"
#include <cmath>
#include <iomanip>

void DomainComputation::initialize(int argc, char *argv[]) {
    
    settings_.loadFromFile(argv[1]);

    partitioning_ = std::make_shared<Partitioning>();
    partitioning_->initialize(settings_.nCells);
    std::array<int,2> nCellsLocal = partitioning_->nCellsLocal();

    // calculate mesh width
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    domain_ = std::make_shared<Domain>(&settings_, partitioning_);
    domain_->readDomainFile(argv[2]);

    // only print for rank 0
    if (partitioning_->ownRankNo() == 0) {

        // pretty print array2d of domain: the right and top faces array and the obstacle array
        Array2D& obstacleMask = *(domain_->obstacleMaskGlobal_);
        Array2D& rightFacesBC = *(domain_->rightFacesBCGlobal_);
        Array2D& topFacesBC = *(domain_->topFacesBCGlobal_);
        int nCellsX = settings_.nCells[0];
        int nCellsY = settings_.nCells[1];
        std::cout << "Obstacle Mask:" << std::endl;
        for (int j = nCellsY - 1; j >= 0; --j) {
            for (int i = 0; i < nCellsX; ++i) {
                std::cout << (obstacleMask(i, j) > 0.5 ? '#' : '.') ;
            }
            std::cout << std::endl;
        }
        std::cout << "Right Faces BC:" << std::endl;
        for (int j = nCellsY - 1; j >= 0; --j) {
            for (int i = 0; i < nCellsX + 1; ++i) {
                double code = rightFacesBC(i, j);
                char marker = domain_->rightFaceMarkerMap().count(code) ? domain_->rightFaceMarkerMap().at(code) : '?';
                std::cout << marker ;
            }
            std::cout << std::endl;
        }
        std::cout << "Top Faces BC:" << std::endl;
        for (int j = nCellsY; j >= 0; --j) {
            for (int i = 0; i < nCellsX; ++i) {
                double code = topFacesBC(i, j);
                char marker = domain_->topFaceMarkerMap().count(code) ? domain_->topFaceMarkerMap().at(code) : '?';
                std::cout << marker ;
            }
            std::cout << std::endl;
        }


        // print the maps of the domain
        std::cout << "Right Face BC Info Map:" << std::endl;
        for (const auto& pair : domain_->rightFaceBCInfoMap()) {
            std::cout << "Code: " << pair.first << ", Info: " << pair.second.toString() << std::endl;
        }
        std::cout << "Top Face BC Info Map:" << std::endl;
        for (const auto& pair : domain_->topFaceBCInfoMap()) {
            std::cout << "Code: " << pair.first << ", Info: " << pair.second.toString() << std::endl;
        }   
        // print char to code maps
        std::cout << "Right Face Marker Map:" << std::endl;
        for (const auto& pair : domain_->rightFaceMarkerMap()) {
            std::cout << "Code: " << pair.first << ", Marker: " << pair.second << std::endl;
        }
        std::cout << "Top Face Marker Map:" << std::endl;
        for (const auto& pair : domain_->topFaceMarkerMap()) {
            std::cout << "Code: " << pair.first << ", Marker: " << pair.second << std::endl;
        }

        // print lists of cells
        std::vector<CellInfo> allCellsInfo = domain_->getInfoListAll();
        std::cout << "All Cells Info List:" << std::endl;
        for (const auto& cellInfo : allCellsInfo) {
            std::cout << cellInfo.toString() << std::endl; 
        }
        std::vector<CellInfo> fluidCellsInfo = domain_->getInfoListFluid();
        std::cout << "Fluid Cells Info List:" << std::endl;
        for (const auto& cellInfo : fluidCellsInfo) {
            std::cout << cellInfo.toString() << std::endl;
        }
        std::vector<CellInfo> redCellsInfo = domain_->getRedListFluid();
        std::cout << "Red Fluid Cells Info List:" << std::endl;
        for (const auto& cellInfo : redCellsInfo) {
            std::cout << cellInfo.toString() << std::endl;
        }

        std::vector<CellInfo> blackCellsInfo = domain_->getBlackListFluid();
        std::cout << "Black Fluid Cells Info List:" << std::endl;
        for (const auto& cellInfo : blackCellsInfo) {
            std::cout << cellInfo.toString() << std::endl;
        }

        std::cout << "Ghost Cells Info List:" << std::endl;
        std::vector<CellInfo> ghostCellsInfo = domain_->getGhostList();
        std::cout << "Ghost Cells Info List length: " << ghostCellsInfo.size() << std::endl;
        for (const auto& cellInfo : ghostCellsInfo) {
            std::cout << cellInfo.toString() << std::endl;
        }

    }




    // create discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(nCellsLocal, meshWidth_, settings_.alpha, partitioning_);
    } else {
        std::cout << "ERROR : Only DonorCell discretization is implemented for domain decomposition!" << std::endl;
    }

    // create pressure solver
    // TODOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
    if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<DomainRBGaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_, domain_);
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
        dt_ = 0.05;

        if (currentTime + dt_ > settings_.endTime - time_eps) {
            dt_ = settings_.endTime - currentTime;
        }


        computePreliminaryVelocities();
        computeRightHandSide();
        computePressure();
        computeVelocities();

        currentTime += dt_;
        iterationCount++;

        printProgress(currentTime, iterationCount);

        // this was the fix!!!
        // if (currentTime >= nOutputs) {
        //     outputWriterParaview_->writeFile(currentTime);
        //     nOutputs = nOutputs + 1;
        // }
        outputWriterParaview_->writeFile(currentTime);
        outputWriterText_->writeFile(currentTime);
    }
//    computeTimeStepWidth();
//    std::cout << "Initial time step width dt: " << dt_ << std::endl;
//    applyInitialBoundaryValues();
//    computePreliminaryVelocities();
//    computeVelocities();
//    communicateGhostCells();
//    computeTimeStepWidth();
//    computePreliminaryVelocities();
//    computeVelocities();

   //print discretization_->u in pretty print
    // if (partitioning_->ownRankNo() == 0) {
    
        // if (partitioning_->ownRankNo() == 0) {
        //     std::cout << "Velocity u field (Rank 0):" << std::endl;
        //     int nCellsX = discretization_->nCells()[0];
        //     int nCellsY = discretization_->nCells()[1];
        //     for (int j=nCellsY+2; j >= 0; --j) {
        //         for (int i = 0; i <= nCellsX + 2; ++i) {
        //             std::cout << std::fixed << std::setprecision(4) << discretization_->u(i, j) << " ";
        //         }
        //         std::cout << std::endl;
        //     }
        // }

        // Print u field for each rank separately with synchronization
//         MPI_Barrier(cartComm_);
//         for (int rank = 0; rank < 4; ++rank) {
//             if (partitioning_->ownRankNo() == rank) {
//                 std::cout << "\nVelocity u field (Rank " << rank << "):" << std::endl;
//                 std::cout << partitioning_->ownPartitionContainsLeftBoundary() << std::endl;
//                 int nCellsX = discretization_->nCells()[0];
//                 int nCellsY = discretization_->nCells()[1];
//                 for (int j=nCellsY+2; j >= 0; --j) {
//                     for (int i = 0; i <= nCellsX + 2; ++i) {
//                         std::cout << std::fixed << std::setprecision(4) << discretization_->u(i, j) << " ";
//                     }
//                     std::cout << std::endl;
//                 }
//                 std::cout.flush();
//             }
//             MPI_Barrier(cartComm_);
//         }
//     // }

//    outputWriterParaview_->writeFile(currentTime);
//    outputWriterText_->writeFile(currentTime);
}

void DomainComputation::printProgress(double &currentTime, int &iterationCount)
{
    if (partitioning_->ownRankNo() == 0)
    {
        if (iterationCount % 10 == 0 || currentTime >= settings_.endTime)
        {
            int percent = static_cast<int>((currentTime / settings_.endTime) * 100);

            // Create progress bar
            const int barWidth = 40;
            int pos = barWidth * currentTime / settings_.endTime;
            std::string progressBar = "[";
            for (int i = 0; i < barWidth; ++i)
            {
                if (i < pos)
                    progressBar += "=";
                else if (i == pos)
                    progressBar += ">";
                else
                    progressBar += " ";
            }
            progressBar += "]";

            std::cout << "\rProgress: " << progressBar << " " << percent << "% | Time: " << currentTime
                      << "/" << settings_.endTime << " | Iter: " << iterationCount << std::flush;
        }
    }
}

void DomainComputation::applyInitialBoundaryValues() {
    //it is sufficient to only go though all cells and only look at the right and top face since then we will go through all faces exactly once. the values set here are only when we have dirichlet BCs directly orthogonally flowing through the face direction. These are onyl set once in the beginning and then never touched again.
    std::vector<CellInfo> fluidCellsInfo = domain_->getInfoListFluid();
    int n = fluidCellsInfo.size();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = fluidCellsInfo[idx];
        if (cellInfo.hasAnyBoundaryFace()) {
            int i = cellInfo.cellIndexPartition[0];
            int j = cellInfo.cellIndexPartition[1];
            // top face
            if (cellInfo.topIsBoundaryFace() && cellInfo.faceTop.dirichletV.has_value()) {
                double faceBCValue = cellInfo.faceTop.dirichletV.value();
                discretization_->v(i, j) = faceBCValue;
                discretization_->g(i, j) = faceBCValue;
                // std::cout << "Setting initial top BC at cell (" << i << ", " << j << ") to " << faceBCValue << std::endl;
            }
            // bottom face
            if (cellInfo.bottomIsBoundaryFace() && cellInfo.faceBottom.dirichletV.has_value()) {
                double faceBCValue = cellInfo.faceBottom.dirichletV.value();
                discretization_->v(i, j - 1) = faceBCValue;
                discretization_->g(i, j - 1) = faceBCValue;
                // std::cout << "Setting initial bottom BC at cell (" << i << ", " << j - 1 << ") to " << faceBCValue << std::endl;
            }
            // right face
            if (cellInfo.rightIsBoundaryFace() && cellInfo.faceRight.dirichletU.has_value()) {
                double faceBCValue = cellInfo.faceRight.dirichletU.value();
                discretization_->u(i, j) = faceBCValue;
                discretization_->f(i, j) = faceBCValue;
                // std::cout << "Setting initial right BC at cell (" << i << ", " << j << ") to " << faceBCValue << std::endl;
            }
            // left face
            if (cellInfo.leftIsBoundaryFace() && cellInfo.faceLeft.dirichletU.has_value()) {
                double faceBCValue = cellInfo.faceLeft.dirichletU.value();
                discretization_->u(i - 1, j) = faceBCValue;
                discretization_->f(i - 1, j) = faceBCValue;
                // std::cout << "Setting initial left BC at cell (" << i - 1 << ", " << j << ") to " << faceBCValue << std::endl;
            }
        }
    }
}


void DomainComputation::communicateGhostCells() {
    const int uIBegin = 1;
    const int uIEnd = discretization_->nCells()[0] + 1;
    const int uJBegin = 1;
    const int uJEnd = discretization_->uJEnd();

    const int vIBegin = 1;
    const int vIEnd = discretization_->vIEnd();
    const int vJBegin = 1;
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
    std::vector<CellInfo> allCellsInfo = domain_->getInfoListFluid();
    std::vector<CellInfo> ghostCellsInfo = domain_->getGhostList();
    allCellsInfo.insert(allCellsInfo.end(), ghostCellsInfo.begin(), ghostCellsInfo.end()); // add ghost cells, but make sure for these only the faces are calculated that are not ghost faces
    int n = allCellsInfo.size();
    double dx = discretization_->dx();
    double dy = discretization_->dy();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = allCellsInfo[idx];
        const int i = cellInfo.cellIndexPartition[0];
        const int j = cellInfo.cellIndexPartition[1];
        // if (cellInfo.hasAnyBoundaryFace()) {
            // check which faces are boundary faces and apply one sided stencils accordingly
        const bool topBC = cellInfo.topIsBoundaryFace();
        const bool rightBC = cellInfo.rightIsBoundaryFace();
        const bool bottomBC = cellInfo.bottomIsBoundaryFace();
        const bool leftBC = cellInfo.leftIsBoundaryFace();

        const double u_i_j = discretization_->u(i,j);
        double u_ip1_j = discretization_->u(i+1,j);
        double u_im1_j = discretization_->u(i-1,j);
        double u_i_jp1 = discretization_->u(i,j+1);
        double u_i_jm1 = discretization_->u(i,j-1);
        double u_im1_jp1 = discretization_->u(i-1,j+1);
        const double v_i_j = discretization_->v(i,j);
        double v_ip1_j = discretization_->v(i+1,j);
        double v_im1_j = discretization_->v(i-1,j);
        double v_i_jm1 = discretization_->v(i,j-1);
        double v_ip1_jm1 = discretization_->v(i+1,j-1);
        double v_i_jp1 = discretization_->v(i,j+1);

        bool calcA = true; // corresponds to f
        bool calcB = true; // correspomds to g

        bool isLeftGhost = i == 1; // these indexes only exist if the cell is a ghost cell
        bool isBottomGhost = j == 1; 
        if (isLeftGhost) {
            // std::cout << "Cell (" << i << ", " << j << ") is a left ghost cell." << std::endl;
            calcB = false; // in left ghost cells we do not calculate g since v there is only set via communication
        }
        if (isBottomGhost) {
            // std::cout << "Cell (" << i << ", " << j << ") is a bottom ghost cell." << std::endl;
            calcA = false; // in bottom ghost cells we do not calculate f since u there is only set via communication
        }

            
        if (rightBC) {
            // std::cout << "Applying right BC at cell (" << i << ", " << j << ")" << std::endl;
            if (cellInfo.faceRight.dirichletU.has_value()) {
                discretization_->f(i,j) = cellInfo.faceRight.dirichletU.value();
                // discretization_->u(i,j) = cellInfo.faceRight.dirichletU.value(); // we have to set the values here because in the parallelcomputation we also first apply the boundaries and then compute the preliminary velocities
                // u_i_j = cellInfo.faceRight.dirichletU.value();
                // TODO SHould we explicitly set u_i_j and its field variable? OR check for assumtions via debuggint that this value is never touched
                calcA = false;
                // std::cout << "  Dirichlet U BC: f(i,j) = " << discretization_->f(i,j) << std::endl;
            } else if (cellInfo.faceRight.neumannU.has_value()) {
                if (cellInfo.faceLeft.neumannU.has_value()) {
                    std::cout << "ERROR: Neumann BCs on both sides of the cell at (" << i << ", " << j << ")" << std::endl;
                    // throw error
                    std::exit(EXIT_FAILURE);
                } else {
                    discretization_->f(i,j) = u_im1_j + cellInfo.faceRight.neumannU.value() * dx;
                    // discretization_->u(i,j) = u_im1_j + cellInfo.faceRight.neumannU.value() * dx;
                    // u_i_j = discretization_->u(i,j);
                    calcA = false;       
                }
                // std::cout << "  Neumann U BC: f(i,j) = " << discretization_->f(i,j) << std::endl;
            }
            if (cellInfo.faceRight.dirichletV.has_value()) {
                v_ip1_j = 2 * cellInfo.faceRight.dirichletV.value() - v_i_j;
                // u_ip1_j does not need to be set since f will not be calculated since if we have a bc to the right, means either dirichlet or neumann for u on the right face, so calcA is false
                // TODO check if this assumtion holds in the domain.cpp
                // std::cout << "  Dirichlet V BC: v_ip1_j = " << v_ip1_j << std::endl;
            }  else if (cellInfo.faceRight.neumannV.has_value()) {
                v_ip1_j = v_i_j + cellInfo.faceRight.neumannV.value() * dx;
                // std::cout << "  Neumann V BC: v_ip1_j = " << v_ip1_j << std::endl;
            }
        }

        if (topBC) {
            if (cellInfo.faceTop.dirichletV.has_value()) {
                discretization_->g(i,j) = cellInfo.faceTop.dirichletV.value();
                calcB = false;
                // std::cout << "  Dirichlet V BC: g(i,j) = " << discretization_->g(i,j) << std::endl;
            } else if (cellInfo.faceTop.neumannV.has_value()) {
                if (cellInfo.faceBottom.neumannV.has_value()) {
                    std::cout << "ERROR: Neumann BCs on both sides of the cell at (" << i << ", " << j << ")" << std::endl;
                    // throw error
                    std::exit(EXIT_FAILURE);
                } else {
                    discretization_->g(i,j) = v_i_jm1 + cellInfo.faceTop.neumannV.value() * dy;
                    calcB = false;
                }
                // std::cout << "  Neumann V BC: g(i,j) = " << discretization_->g(i,j) << std::endl;
            }
            // std::cout << "Applying top BC at cell (" << i << ", " << j << ")" << std::endl;
            if (cellInfo.faceTop.dirichletU.has_value()) {
                u_i_jp1 = 2.0 * cellInfo.faceTop.dirichletU.value() - u_i_j;
                // by settng u_i_j in the right face bc above first when we do dirichlet there, we effectively prioritize right  faces over top faces if the conditions conflict.
                // because if we have neumann u at the top and dirichlet u at the rigth it doesnt work so we first set the right face u value and then use it here to set the top ghost value
                // std::cout << "  Dirichlet U BC: u_i_jp1 = " << u_i_jp1 << std::endl;
            }  else if (cellInfo.faceTop.neumannU.has_value()) {
                u_i_jp1 = u_i_j + cellInfo.faceTop.neumannU.value() * dy;
                // std::cout << "  Neumann U BC: u_i_jp1 = " << u_i_jp1 << std::endl;
            }
        }
            
        if (leftBC) {
            // std::cout << "Applying left BC at cell (" << i << ", " << j << ")" << std::endl;
            if (cellInfo.faceLeft.dirichletU.has_value()) {
                // left boundary has Dirichlet for u, ghost value already set in applyInitialBoundaryValues
                // u_im1_j = cellInfo.faceLeft.dirichletU.value(); //??????
                // discretization_->f(i-1,j) = u_im1_j;
                // continue;
            } else if (cellInfo.faceLeft.neumannU.has_value()) { 
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
            // std::cout << "Applying bottom BC at cell (" << i << ", " << j << ")" << std::endl;
            if (cellInfo.faceBottom.dirichletV.has_value()) {
                // bottom boundary has Dirichlet for v, ghost value already set in applyInitialBoundaryValues
                // v_i_jm1 = 2.0 * cellInfo.faceBottom.dirichletV.value() - v_i_j; //?????
                // continue;
            } else if (cellInfo.faceBottom.neumannV.has_value()) {
                v_i_jm1 = v_i_j + cellInfo.faceBottom.neumannV.value() * dy;
                discretization_->g(i,j-1) = v_i_jm1; // set g to the neumann value
            }
            if (cellInfo.faceBottom.dirichletU.has_value()) {
                u_i_jm1 = 2.0 * cellInfo.faceBottom.dirichletU.value() - u_i_j;
            }  else if (cellInfo.faceBottom.neumannU.has_value()) {
                u_i_jm1 = u_i_j + cellInfo.faceBottom.neumannU.value() * dy;     
            }
            
        }
            // if (i == 2 && j == 11 && partitioning_->ownRankNo() == 1) {
            //     std::cout << calcA << "BBBBBBBBBBBBBBBBBBBBBBBBBis calcA for cell (2,11) in rank 1" << std::endl;
                
            // }

            
        if (calcA) {
            // if (topBC) {
            //     std::cout << "  Using one-sided stencil for D2uDy2 at top BC" << std::endl;
            //     std::cout << "    u_i_jp1 = " << u_i_jp1 << std::endl;
            //     std::cout << calcA << std::endl;
            // }
            double A_ij = 1 / settings_.re * (computeD2uDx2(u_ip1_j, u_i_j, u_im1_j) + computeD2uDy2(u_i_jp1, u_i_j, u_i_jm1)) - computeDu2Dx(u_i_j, u_im1_j, u_ip1_j) - computeDuvDy(u_i_j, u_i_jp1, u_i_jm1, v_i_j, v_ip1_j, v_i_jm1, v_ip1_jm1) + settings_.g[0];
            discretization_->f(i,j) = u_i_j + A_ij * dt_;
            // std::cout << "Applying A_ij calculation at cell (" << i << ", " << j << ") in rank " << partitioning_->ownRankNo() << std::endl;
            // std::cout << "  Calculated f(i,j) = " << discretization_->f(i,j) << "with A_ij = " << A_ij << "and dt = " << dt_ << std::endl;
        }

        if (calcB) {
            double B_ij = 1 / settings_.re * (computeD2vDx2(v_ip1_j, v_i_j, v_im1_j) + computeD2vDy2(v_i_jp1, v_i_j, v_i_jm1)) - computeDuvDx(u_i_j, u_i_jp1, u_im1_j, u_im1_jp1, v_i_j, v_ip1_j, v_im1_j) - computeDv2Dy(v_i_j, v_i_jp1, v_i_jm1) + settings_.g[1];
            discretization_->g(i,j) = v_i_j + B_ij * dt_;   
        }      

        // // } else {
        //     // inner cell, normal stencil
            // double A_ij = 1 / settings_.re * ( discretization_->computeD2uDx2(i,j) + discretization_->computeD2uDy2(i,j)) - discretization_->computeDu2Dx(i,j) - discretization_->computeDuvDy(i,j) + settings_.g[0];
        //     discretization_->f(i,j) = discretization_->u(i,j) + A_ij * dt_;

            // double B_ij = 1 / settings_.re * ( discretization_->computeD2vDx2(i,j) + discretization_->computeD2vDy2(i,j)) - discretization_->computeDuvDx(i,j) - discretization_->computeDv2Dy(i,j) + settings_.g[1];
        //     discretization_->g(i,j) = discretization_->v(i,j) + B_ij * dt_;
        //     // if (i == 2 && j == 11 && partitioning_->ownRankNo() == 1) {
        //     //     std::cout << "AAAAAAAAAAAAAAAAAAAAAAAAAAAAis calcA for cell (2,11) in rank 1" << std::endl;
                
        //     // }
        // // }
        

        // if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        //     if (i == 2 && cellInfo.leftIsBoundaryFace() == false) {
        //         std::cout << "Computing left ghost cell for cell (" << i << ", " << j << ")" << std::endl;
        //         // compute u at left ghost cell if left boundary is in this partition
        //         // double A_i_minus_1_j = 1 / settings_.re * ( discretization_->computeD2uDx2(i-1,j) + discretization_->computeD2uDy2(i-1,j)) - discretization_->computeDu2Dx(i-1,j) - discretization_->computeDuvDy(i-1,j) + settings_.g[0];
        //         // discretization_->f(i-1,j) = discretization_->u(i-1,j) + A_i_minus_1_j * dt_;
        //         // everything same as above but with i-1
        //         double u_im2_j = discretization_->u(i-2,j);
        //         double u_im1_jm1 = discretization_->u(i-1,j-1);
        //         double v_im1_jm1 = discretization_->v(i-1,j-1);
        //         double A_im1j = 1 / settings_.re * (computeD2uDx2(u_i_j, u_im1_j, u_im2_j) + computeD2uDy2(u_im1_jp1, u_im1_j, u_im1_jm1)) - computeDu2Dx(u_im1_j, u_im2_j, u_i_j) - computeDuvDy(u_im1_j, u_im1_jp1, u_im1_jm1, v_im1_j, v_i_j, v_im1_jm1, v_i_jm1) + settings_.g[0];
        //         std::cout << "  A_im1j = " << A_im1j << "for left ghost cell at (" << i-1 << ", " << j << ") in rank " << partitioning_->ownRankNo() << std::endl;
        //         // Print debug info for left ghost cell u calculation only for rank 3 and specific cell
        //         if (partitioning_->ownRankNo() == 3 && i == 2 && j == 11) {
        //             std::cout << "\n=== Left Ghost Cell U Calculation Debug (i=" << i << ", j=" << j << ", Rank " << partitioning_->ownRankNo() << ") ===" << std::endl;
        //             std::cout << "u values:" << std::endl;
        //             std::cout << "  u_im2_j = u(" << i-2 << "," << j << ") = " << u_im2_j << std::endl;
        //             std::cout << "  u_im1_j = u(" << i-1 << "," << j << ") = " << u_im1_j << std::endl;
        //             std::cout << "  u_i_j = u(" << i << "," << j << ") = " << u_i_j << std::endl;
        //             std::cout << "  u_im1_jm1 = u(" << i-1 << "," << j-1 << ") = " << u_im1_jm1 << std::endl;
        //             std::cout << "  u_im1_jp1 = u(" << i-1 << "," << j+1 << ") = " << u_im1_jp1 << std::endl;
        //             std::cout << "v values:" << std::endl;
        //             std::cout << "  v_im1_j = v(" << i-1 << "," << j << ") = " << v_im1_j << std::endl;
        //             std::cout << "  v_i_j = v(" << i << "," << j << ") = " << v_i_j << std::endl;
        //             std::cout << "  v_im1_jm1 = v(" << i-1 << "," << j-1 << ") = " << v_im1_jm1 << std::endl;
        //             std::cout << "  v_i_jm1 = v(" << i << "," << j-1 << ") = " << v_i_jm1 << std::endl;
        //             double d2u_dx2 = computeD2uDx2(u_i_j, u_im1_j, u_im2_j);
        //             double d2u_dy2 = computeD2uDy2(u_im1_jp1, u_im1_j, u_im1_jm1);
        //             double du2_dx = computeDu2Dx(u_im1_j, u_im2_j, u_i_j);
        //             double duv_dy = computeDuvDy(u_im1_j, u_im1_jp1, u_im1_jm1, v_im1_j, v_i_j, v_im1_jm1, v_i_jm1);
        //             std::cout << "Stencil values:" << std::endl;
        //             std::cout << "  d2u_dx2 = " << d2u_dx2 << std::endl;
        //             std::cout << "  d2u_dy2 = " << d2u_dy2 << std::endl;
        //             std::cout << "  du2_dx = " << du2_dx << std::endl;
        //             std::cout << "  duv_dy = " << duv_dy << std::endl;
        //             std::cout << "Parameters:" << std::endl;
        //             std::cout << "  settings_.re = " << settings_.re << std::endl;
        //             std::cout << "  settings_.g[0] = " << settings_.g[0] << std::endl;
        //             std::cout << "  dt_ = " << dt_ << std::endl;
        //         }
        //         discretization_->f(i-1,j) = u_im1_j + A_im1j * dt_;
        //     }
        // }
        // if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        //     if (j == 2 && cellInfo.bottomIsBoundaryFace() == false) {
        //         // std::cout << "Computing bottom ghost cell for cell (" << i << ", " << j << ")" << std::endl;
        //         // compute v at bottom ghost cell if bottom boundary is in this partition
        //         double v_i_jm2 = discretization_->v(i,j-2);
        //         double v_im1_jm1 = discretization_->v(i-1,j-1);
        //         double u_im1_jm1 = discretization_->u(i-1,j-1);
        //         double B_ijm1 = 1 / settings_.re * (computeD2vDx2(v_ip1_jm1, v_i_jm1, v_im1_jm1) + computeD2vDy2(v_i_j, v_i_jm1, v_i_jm2)) - computeDuvDx(u_i_jm1, u_i_j, u_im1_jm1, u_im1_j, v_i_jm1, v_ip1_jm1, v_im1_jm1) - computeDv2Dy(v_i_jm1, v_i_j, v_i_jm2) + settings_.g[1];
        //         discretization_->g(i,j-1) = v_i_j + B_ijm1 * dt_;   
        //     }
        // }
    }
}

void DomainComputation::computeRightHandSide() {
    std::vector<CellInfo> allCellsInfo = domain_->getInfoListFluid();
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
    std::vector<CellInfo> fluidCellsInfo = domain_->getInfoListFluid();
    int n = fluidCellsInfo.size();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = fluidCellsInfo[idx];
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];

        double p_i_j = discretization_->p(i,j);
        double p_ip1_j = discretization_->p(i+1,j);
        double p_i_jp1 = discretization_->p(i,j+1);


        if (!cellInfo.rightIsBoundaryFace()){
            // std::cout << "Updating u velocity at cell (" << i << ", " << j << ") with f: " << discretization_->f(i,j) << " and pressure gradient: " << computeDpDx(p_ip1_j, p_i_j) << std::endl;
            discretization_->u(i,j) = discretization_->f(i,j) - dt_ * computeDpDx(p_ip1_j, p_i_j);
            // std::cout << "  New u(i,j) = " << discretization_->u(i,j) << std::endl;
        } else if (cellInfo.faceRight.neumannU.has_value()) {
            // right boundary face
            // std::cout << "Applying Neumann BC for u velocity at right boundary face of cell (" << i << ", " << j << ") with f: " << discretization_->f(i,j) << std::endl;
            discretization_->u(i,j) = discretization_->f(i,j);
        }

        if (!cellInfo.topIsBoundaryFace()){
            discretization_->v(i,j) = discretization_->g(i,j) - dt_ * computeDpDy(p_i_jp1, p_i_j);
        } else if (cellInfo.faceTop.neumannV.has_value()) {
            // top boundary face
            discretization_->v(i,j) = discretization_->g(i,j);
        }

        if (!partitioning_->ownPartitionContainsLeftBoundary()) {
            if (!cellInfo.leftIsBoundaryFace() && i-1 == discretization_->uIBegin()) {
                // std::cout << "Computing left ghost u velocity for cell (" << i << ", " << j << ") with f: " << discretization_->f(i-1,j) << " and pressure gradient: " << computeDpDx(p_i_j, discretization_->p(i-1,j)) << std::endl;
                discretization_->u(i-1,j) = discretization_->f(i-1,j) - dt_ * computeDpDx(p_i_j, discretization_->p(i-1,j));
            } else if (cellInfo.faceLeft.neumannU.has_value()) {
                // left boundary face
                // std::cout << "Applying Neumann BC for u velocity at left boundary face of cell (" << i << ", " << j << ") with f: " << discretization_->f(i-1,j) << std::endl;
                discretization_->u(i-1,j) = discretization_->f(i-1,j);
            }
        }

        if (!partitioning_->ownPartitionContainsBottomBoundary()) {
            if (!cellInfo.bottomIsBoundaryFace() && j-1 == discretization_->vJBegin()) {
                discretization_->v(i,j-1) = discretization_->g(i,j-1) - dt_ * computeDpDy(p_i_j, discretization_->p(i,j-1));
            } else if (cellInfo.faceBottom.neumannV.has_value()) {
                // bottom boundary face
                discretization_->v(i,j-1) = discretization_->g(i,j-1);
            }
        }

        if (partitioning_->ownPartitionContainsTopBoundary()) {
            if (j == discretization_->nCells()[1] + 1) {
                // top has to be a ghost cell. we have to set u there since paraview output writer needs that ghost u value to interpolate u on the boundary
                if (cellInfo.faceTop.dirichletU.has_value()) {
                    double faceBCValue = cellInfo.faceTop.dirichletU.value();
                    double u_i_j = discretization_->u(i,j);
                    discretization_->u(i,j+1) = 2.0 * faceBCValue - u_i_j;
                } 
            }
        }
        if (partitioning_->ownPartitionContainsLeftBoundary()) {
            if (i == 2) {
                // left has to be a ghost cell. we have to set v there since paraview output writer needs that ghost v value to interpolate v on the boundary
                if (cellInfo.faceLeft.dirichletV.has_value()) {
                    double faceBCValue = cellInfo.faceLeft.dirichletV.value();
                    double v_i_j = discretization_->v(i,j);
                    discretization_->v(i-1,j) = 2.0 * faceBCValue - v_i_j;
                } 
            }   
        }
        if (partitioning_->ownPartitionContainsBottomBoundary()) {
            if (j == 2) {
                // bottom has to be a ghost cell. we have to set u there since paraview output writer needs that ghost u value to interpolate u on the boundary
                if (cellInfo.faceBottom.dirichletU.has_value()) {
                    double faceBCValue = cellInfo.faceBottom.dirichletU.value();
                    double u_i_j = discretization_->u(i,j);
                    discretization_->u(i,j-1) = 2.0 * faceBCValue - u_i_j;
                } 
                // ???? MÜSSEN WIR NEUMANN AUCH BEACHTEN? TODO
            }   
        }
        if (partitioning_->ownPartitionContainsRightBoundary()) {
            if (i == discretization_->nCells()[0] + 1) {
                // right has to be a ghost cell. we have to set v there since paraview output writer needs that ghost v value to interpolate v on the boundary
                if (cellInfo.faceRight.dirichletV.has_value()) {
                    double faceBCValue = cellInfo.faceRight.dirichletV.value();
                    double v_i_j = discretization_->v(i,j);
                    discretization_->v(i+1,j) = 2.0 * faceBCValue - v_i_j;
                } 
            }
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