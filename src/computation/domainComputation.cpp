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
        obstacleMask.prettyPrintArray2D();
        std::cout << "Right Faces BC:" << std::endl;
        rightFacesBC.prettyPrintArray2D();
        std::cout << "Top Faces BC:" << std::endl;
        topFacesBC.prettyPrintArray2D();


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
    // TODO
    if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<DomainRBGaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_, domain_);
    } else {
        std::cerr << "Error: Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(EXIT_FAILURE);
    }

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
    
    communicateGhostCells();

    while (currentTime < settings_.endTime - time_eps) {

        computeTimeStepWidth();

        if (currentTime + dt_ > settings_.endTime - time_eps) {
            dt_ = settings_.endTime - currentTime;
        }
        if (partitioning_->ownRankNo() == 0) {
            std::cout << dt_ << std::endl;
        }
        


        computePreliminaryVelocities();
        computeRightHandSide();

        MPI_Barrier(cartComm_);
        if (partitioning_->ownRankNo() == 3 && iterationCount < 2) {
            std::cout << "u before pressure solve:" << std::endl;
            discretization_->u().printAsArray();
            std::cout << "f before pressure solve:" << std::endl;
            discretization_->f().printAsArray();
            std::cout << "p before pressure solve:" << std::endl;
            discretization_->p().printAsArray();
        }
        MPI_Barrier(cartComm_);

        computePressure();

        MPI_Barrier(cartComm_);
        if (partitioning_->ownRankNo() == 3 && iterationCount < 2) {
            std::cout << "p after pressure solve:" << std::endl;
            discretization_->p().printAsArray();
        }
        MPI_Barrier(cartComm_);

        computeVelocities();

        MPI_Barrier(cartComm_);
        if (partitioning_->ownRankNo() == 3 && iterationCount < 2) {
            std::cout << "u after calculating u:" << std::endl;
            discretization_->u().printAsArray();
        }
        MPI_Barrier(cartComm_);

        currentTime += dt_;
        iterationCount++;

        printProgress(currentTime, iterationCount);

        // this was the fix!!!
        // if (currentTime >= nOutputs) {
        //     outputWriterParaview_->writeFile(currentTime);
        //     nOutputs = nOutputs + 1;
        // }
        communicateGhostCells();

        outputWriterParaview_->writeFile(currentTime);
        outputWriterText_->writeFile(currentTime);
    }
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
    std::vector<CellInfo> fluidCellsInfo = domain_->getInfoListFluid();
    int n = fluidCellsInfo.size();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = fluidCellsInfo[idx];
        if (cellInfo.hasAnyBoundaryFace()) {
            const int i = cellInfo.cellIndexPartition[0];
            const int j = cellInfo.cellIndexPartition[1];

            if (cellInfo.faceTop.isBoundaryFace() && !cellInfo.faceRight.isBoundaryFace() && !cellInfo.faceLeft.isBoundaryFace()) {
                if (cellInfo.faceTop.dirichletV.has_value()) {
                    double faceBCValue = cellInfo.faceTop.dirichletV.value();
                    discretization_->v(i, j) = faceBCValue;
                    discretization_->g(i, j) = faceBCValue;
                }
                if (cellInfo.faceTop.dirichletU.has_value()) {
                    discretization_->u(i, j+1) = 2*cellInfo.faceTop.dirichletU.value() - discretization_->u(i, j); // mirror value for u at top face
                    // discretization_->u(i-1, j+1) = 2*cellInfo.faceTop.dirichletU.value() - discretization_->u(i, j); // mirror value for u at top face
                }
            }
            
            if (cellInfo.faceBottom.isBoundaryFace() && !cellInfo.faceRight.isBoundaryFace() && !cellInfo.faceLeft.isBoundaryFace()) {
                if (cellInfo.faceBottom.dirichletV.has_value()) {
                    double faceBCValue = cellInfo.faceBottom.dirichletV.value();
                    discretization_->v(i, j - 1) = faceBCValue;
                    discretization_->g(i, j - 1) = faceBCValue;
                }
                if (cellInfo.faceBottom.dirichletU.has_value()) {
                    discretization_->u(i, j-1) = 2*cellInfo.faceBottom.dirichletU.value() - discretization_->u(i, j); // mirror value for u at bottom face
                }
            }
            
            if (cellInfo.faceRight.isBoundaryFace()) {
                if (cellInfo.faceRight.dirichletU.has_value()) {
                    double faceBCValue = cellInfo.faceRight.dirichletU.value();
                    discretization_->u(i, j) = faceBCValue;
                    discretization_->f(i, j) = faceBCValue;
                }
                if (cellInfo.faceRight.dirichletV.has_value()) {
                    discretization_->v(i + 1, j) = 2*cellInfo.faceRight.dirichletV.value() - discretization_->v(i, j); // mirror value for v at right face
                }
            }
           
            if (cellInfo.faceLeft.isBoundaryFace()) {
                if (cellInfo.faceLeft.dirichletU.has_value()) {
                    double faceBCValue = cellInfo.faceLeft.dirichletU.value();
                    discretization_->u(i - 1, j) = faceBCValue;
                    discretization_->f(i - 1, j) = faceBCValue;
                }
                if (cellInfo.faceLeft.dirichletV.has_value()) {
                    discretization_->v(i - 1, j) = 2*cellInfo.faceLeft.dirichletV.value() - discretization_->v(i, j); // mirror value for v at left face
                }
            }
        }
    }
}


void DomainComputation::communicateGhostCells() {
    const int iBegin = 1;
    const int iEnd = discretization_->nCells()[0] + 2;
    const int jBegin = 1;
    const int jEnd = discretization_->nCells()[1] + 2;
    const int l = iEnd - iBegin + 1;

    // buffers for sending and receiving data
    std::vector<double> sendBufferTopU(l, 0.0);
    std::vector<double> sendBufferTopV(l, 0.0);
    std::vector<double> sendBufferBottomU(l, 0.0);
    std::vector<double> sendBufferBottomV(l, 0.0);
    std::vector<double> sendBufferLeftU(l, 0.0);
    std::vector<double> sendBufferLeftV(l, 0.0);
    std::vector<double> sendBufferRightU(l, 0.0);
    std::vector<double> sendBufferRightV(l, 0.0);

    // buffors for receiving data
    std::vector<double> recvBufferTopU(l, 0.0);
    std::vector<double> recvBufferTopV(l, 0.0);
    std::vector<double> recvBufferBottomU(l, 0.0);
    std::vector<double> recvBufferBottomV(l, 0.0);
    std::vector<double> recvBufferLeftU(l, 0.0);
    std::vector<double> recvBufferLeftV(l, 0.0);
    std::vector<double> recvBufferRightU(l, 0.0);
    std::vector<double> recvBufferRightV(l, 0.0);

    MPI_Request requestsSendTopU, requestsSendTopV, requestsRecvTopU, requestsRecvTopV;
    MPI_Request requestsSendBottomU, requestsSendBottomV, requestsRecvBottomU, requestsRecvBottomV;
    MPI_Request requestsSendLeftU, requestsSendLeftV, requestsRecvLeftU, requestsRecvLeftV;
    MPI_Request requestsSendRightU, requestsSendRightV, requestsRecvRightU, requestsRecvRightV;
    const int TAG_U = 0;
    const int TAG_V = 1;

    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        // otherwise communicate with the top neighbour
        for (int i = iBegin; i <= iEnd; i++) {
            sendBufferTopU[i - iBegin] = discretization_->u(i,jEnd - 1);
            sendBufferTopV[i - iBegin] = discretization_->v(i,jEnd - 2);
        }
        // instantiate non-blocking sends and receives
        MPI_Isend(sendBufferTopU.data(), sendBufferTopU.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), TAG_U, cartComm_, &requestsSendTopU);
        MPI_Isend(sendBufferTopV.data(), sendBufferTopV.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), TAG_V, cartComm_, &requestsSendTopV);

        MPI_Irecv(recvBufferTopU.data(), recvBufferTopU.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), TAG_U, cartComm_, &requestsRecvTopU);
        MPI_Irecv(recvBufferTopV.data(), recvBufferTopV.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), TAG_V, cartComm_, &requestsRecvTopV);
    }
    
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = iBegin; i <= iEnd; i++) {
            sendBufferBottomU[i - iBegin] = discretization_->u(i,jBegin + 1);
            sendBufferBottomV[i - iBegin] = discretization_->v(i,jBegin + 1); // +1 because we have two layers of gjost cells at the bottom (just like at the left)
        }
        MPI_Isend(sendBufferBottomU.data(), sendBufferBottomU.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_U, cartComm_, &requestsSendBottomU);
        MPI_Isend(sendBufferBottomV.data(), sendBufferBottomV.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_V, cartComm_, &requestsSendBottomV);

        MPI_Irecv(recvBufferBottomU.data(), recvBufferBottomU.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_U, cartComm_, &requestsRecvBottomU);
        MPI_Irecv(recvBufferBottomV.data(), recvBufferBottomV.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_V, cartComm_, &requestsRecvBottomV);
    }

    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = jBegin; j <= jEnd; j++) {
            sendBufferLeftU[j - jBegin] = discretization_->u(iBegin + 1,j); // +1 because we have two layers of ghost cells at the left (just like at the bottom)
            sendBufferLeftV[j - jBegin] = discretization_->v(iBegin + 1,j);
        }
        MPI_Isend(sendBufferLeftU.data(), sendBufferLeftU.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), TAG_U, cartComm_, &requestsSendLeftU);
        MPI_Isend(sendBufferLeftV.data(), sendBufferLeftV.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), TAG_V, cartComm_, &requestsSendLeftV);
        
        MPI_Irecv(recvBufferLeftU.data(), recvBufferLeftU.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), TAG_U, cartComm_, &requestsRecvLeftU);
        MPI_Irecv(recvBufferLeftV.data(), recvBufferLeftV.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), TAG_V, cartComm_, &requestsRecvLeftV);
    }

    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = jBegin; j <= jEnd; j++) {
            sendBufferRightU[j - jBegin] = discretization_->u(iEnd - 2,j);
            sendBufferRightV[j - jBegin] = discretization_->v(iEnd - 1,j);
        }
        MPI_Isend(sendBufferRightU.data(), sendBufferRightU.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), TAG_U, cartComm_, &requestsSendRightU);
        MPI_Isend(sendBufferRightV.data(), sendBufferRightV.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), TAG_V, cartComm_, &requestsSendRightV);

        MPI_Irecv(recvBufferRightU.data(), recvBufferRightU.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), TAG_U, cartComm_, &requestsRecvRightU);
        MPI_Irecv(recvBufferRightV.data(), recvBufferRightV.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), TAG_V, cartComm_, &requestsRecvRightV);
    }

    // wait for all communications to finish and set ghost values
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestsRecvTopU, MPI_STATUS_IGNORE);
        MPI_Wait(&requestsRecvTopV, MPI_STATUS_IGNORE);
        for (int i = iBegin; i <= iEnd; i++) {
            discretization_->u(i,jEnd) = recvBufferTopU[i - iBegin];
            discretization_->v(i,jEnd) = recvBufferTopV[i - iBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestsRecvBottomU, MPI_STATUS_IGNORE);
        MPI_Wait(&requestsRecvBottomV, MPI_STATUS_IGNORE);
        for (int i = iBegin; i <= iEnd; i++) {
            discretization_->u(i,jBegin) = recvBufferBottomU[i - iBegin];
            discretization_->v(i,jBegin - 1) = recvBufferBottomV[i - iBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestsRecvLeftU, MPI_STATUS_IGNORE);
        MPI_Wait(&requestsRecvLeftV, MPI_STATUS_IGNORE);
        for (int j = jBegin; j <= jEnd; j++) {
            discretization_->u(iBegin - 1,j) = recvBufferLeftU[j - jBegin];
            discretization_->v(iBegin,j) = recvBufferLeftV[j - jBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestsRecvRightU, MPI_STATUS_IGNORE);
        MPI_Wait(&requestsRecvRightV, MPI_STATUS_IGNORE);
        for (int j = jBegin; j <= jEnd; j++) {
            discretization_->u(iEnd,j) = recvBufferRightU[j - jBegin];
            discretization_->v(iEnd,j) = recvBufferRightV[j - jBegin];
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

        bool calcA = !(j == 1); // corresponds to f; is false if the cell is a ghost cell below the partition bc then only g is calculated
        bool calcB = !(i == 1); // correspomds to g; is false if the cell is a ghost cell left of the partition bc then only f is calculated

            
        if (cellInfo.faceRight.isBoundaryFace()) {
            if (cellInfo.faceRight.dirichletU.has_value()) {
                discretization_->f(i,j) = cellInfo.faceRight.dirichletU.value();
                calcA = false;
            } else if (cellInfo.faceRight.neumannU.has_value()) {
                if (cellInfo.faceLeft.neumannU.has_value()) {
                    std::cout << "ERROR: Neumann BCs on both sides of the cell at (" << i << ", " << j << ")" << std::endl;
                    std::exit(EXIT_FAILURE);
                } else {
                    discretization_->f(i,j) = u_im1_j + cellInfo.faceRight.neumannU.value() * dx;
                    calcA = false;       
                }
            }
            if (cellInfo.faceRight.dirichletV.has_value()) {
                v_ip1_j = 2 * cellInfo.faceRight.dirichletV.value() - v_i_j;
            }  else if (cellInfo.faceRight.neumannV.has_value()) {
                v_ip1_j = v_i_j + cellInfo.faceRight.neumannV.value() * dx;
            }
        }

        if (cellInfo.faceTop.isBoundaryFace()) {
            if (cellInfo.faceTop.dirichletV.has_value()) {
                discretization_->g(i,j) = cellInfo.faceTop.dirichletV.value();
                calcB = false;
            } else if (cellInfo.faceTop.neumannV.has_value()) {
                if (cellInfo.faceBottom.neumannV.has_value()) {
                    std::cout << "ERROR: Neumann BCs on both sides of the cell at (" << i << ", " << j << ")" << std::endl;
                    std::exit(EXIT_FAILURE);
                } else {
                    discretization_->g(i,j) = v_i_jm1 + cellInfo.faceTop.neumannV.value() * dy;
                    calcB = false;
                }
            }
            if (cellInfo.faceTop.dirichletU.has_value()) {
                u_i_jp1 = 2.0 * cellInfo.faceTop.dirichletU.value() - u_i_j;
            }  else if (cellInfo.faceTop.neumannU.has_value()) {
                u_i_jp1 = u_i_j + cellInfo.faceTop.neumannU.value() * dy;
            }
        }
            
        if (cellInfo.faceLeft.isBoundaryFace()) {
            if (cellInfo.faceLeft.dirichletU.has_value()) {
                // we trust that applyInitialBoundaryValues has already set the ghost value for dirichlet u at left face
                if (discretization_->u(i-1,j) != cellInfo.faceLeft.dirichletU.value()) {
                    std::cout << "ERROR: Dirichlet u BC at left face of cell (" << i << ", " << j << ") not properly set in applyInitialBoundaryValues!" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else if (cellInfo.faceLeft.neumannU.has_value()) { 
                u_im1_j = u_i_j + cellInfo.faceLeft.neumannU.value() * dx;
                discretization_->f(i-1,j) = u_im1_j; // set f to the neumann value (corresponds to a solid cell bc neumann left means solid obstacle to the left)
            } 
            if (cellInfo.faceLeft.dirichletV.has_value()) {
                v_im1_j = 2.0 * cellInfo.faceLeft.dirichletV.value() - v_i_j;
                if (partitioning_->ownPartitionContainsLeftBoundary() && (i == 1)) {
                    discretization_->g(i-1,j) = v_im1_j; // set g to the dirichlet value (corresponds to a solid cell bc dirichlet left means solid obstacle to the left)
                }
            }  else if (cellInfo.faceLeft.neumannV.has_value()) {
                v_im1_j = v_i_j + cellInfo.faceLeft.neumannV.value() * dx;
                if (partitioning_->ownPartitionContainsLeftBoundary() && (i == 1)) {
                    discretization_->g(i-1,j) = v_im1_j; // set g to the neumann value (corresponds to a solid cell bc neumann left means solid obstacle to the left)
                }
            }
        }
        if (cellInfo.faceBottom.isBoundaryFace()) {
            if (cellInfo.faceBottom.dirichletV.has_value()) {
                // we trust that applyInitialBoundaryValues has already set the ghost value for dirichlet v at bottom face
                if (discretization_->v(i,j-1) != cellInfo.faceBottom.dirichletV.value()) {
                    std::cout << "ERROR: Dirichlet v BC at bottom face of cell (" << i << ", " << j << ") not properly set in applyInitialBoundaryValues!" << std::endl;
                    std::exit(EXIT_FAILURE);
                }
            } else if (cellInfo.faceBottom.neumannV.has_value()) {
                v_i_jm1 = v_i_j + cellInfo.faceBottom.neumannV.value() * dy;
                discretization_->g(i,j-1) = v_i_jm1; // set g to the neumann value
            }
            if (cellInfo.faceBottom.dirichletU.has_value()) {
                u_i_jm1 = 2.0 * cellInfo.faceBottom.dirichletU.value() - u_i_j;
                if (partitioning_->ownPartitionContainsBottomBoundary() && (j == 1)) {
                    discretization_->f(i,j-1) = u_i_jm1; // set f to the dirichlet value (corresponds to a solid cell bc dirichlet bottom means solid obstacle to the bottom)
                }
            }  else if (cellInfo.faceBottom.neumannU.has_value()) {
                u_i_jm1 = u_i_j + cellInfo.faceBottom.neumannU.value() * dy;     
                if (partitioning_->ownPartitionContainsBottomBoundary() && (j == 1)) {
                    discretization_->f(i,j-1) = u_i_jm1; // set f to the neumann value (corresponds to a solid cell bc neumann bottom means solid obstacle to the bottom)
                }
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
    std::vector<CellInfo> ghostCellsInfo = domain_->getGhostList();
    fluidCellsInfo.insert(fluidCellsInfo.end(), ghostCellsInfo.begin(), ghostCellsInfo.end()); // add ghost cells, but make sure for these only the faces are calculated that are not ghost faces
    int n = fluidCellsInfo.size();
    for (int idx = 0; idx < n; idx++) {
        CellInfo cellInfo = fluidCellsInfo[idx];
        int i = cellInfo.cellIndexPartition[0];
        int j = cellInfo.cellIndexPartition[1];

        double p_i_j = discretization_->p(i,j);
        double p_ip1_j = discretization_->p(i+1,j);
        double p_i_jp1 = discretization_->p(i,j+1);

        bool isLeftGhost = i == 1; // these indexes only exist if the cell is a ghost cell
        bool isBottomGhost = j == 1; 

        if (!isBottomGhost) {
            // as long as the cell is not a ghost cell below the partition boundary we can compute u
            if (!cellInfo.faceRight.isBoundaryFace()){
                // if it has no boundary face to the right, then p_ip1_j is valid
                discretization_->u(i,j) = discretization_->f(i,j) - dt_ * computeDpDx(p_ip1_j, p_i_j);
            } else if (cellInfo.faceRight.neumannU.has_value()) {
                // right boundary face: p_ip1_j would be the same as p_i_j, so pressure gradient is zero
                discretization_->u(i,j) = discretization_->f(i,j);
            }
        }
        

        if (!isLeftGhost) {
            // as long as the cell is not a ghost cell left of the partition boundary we can compute v 
            if (!cellInfo.faceTop.isBoundaryFace()){
                // if it has no boundary face at the top, then p_i_jp1 is valid
                discretization_->v(i,j) = discretization_->g(i,j) - dt_ * computeDpDy(p_i_jp1, p_i_j);
            } else if (cellInfo.faceTop.neumannV.has_value()) {
                // top boundary face: pressure gradient is zero
                discretization_->v(i,j) = discretization_->g(i,j);
            }
        }

        // TODO check for sanity
        if (partitioning_->ownPartitionContainsTopBoundary() && !cellInfo.faceRight.isBoundaryFace()) {
            if (j == discretization_->nCells()[1] + 1) {
                // top has to be a ghost cell. we have to set u there since paraview output writer needs that ghost u value to interpolate u on the boundary
                if (cellInfo.faceTop.dirichletU.has_value() ) {
                    double faceBCValue = cellInfo.faceTop.dirichletU.value();
                    double u_i_j = discretization_->u(i,j);
                    discretization_->u(i,j+1) = 2.0 * faceBCValue - u_i_j;
                } else if (cellInfo.faceTop.neumannU.has_value()) {
                    // set ghost value based on neumann condition
                    double neumannValue = cellInfo.faceTop.neumannU.value();
                    double u_i_j = discretization_->u(i,j);
                    discretization_->u(i,j+1) = u_i_j + neumannValue * discretization_->dy();
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
                } else if (cellInfo.faceLeft.neumannV.has_value()) {
                    // set ghost value based on neumann condition
                    double neumannValue = cellInfo.faceLeft.neumannV.value();
                    double v_i_j = discretization_->v(i,j);
                    discretization_->v(i-1,j) = v_i_j + neumannValue * discretization_->dx();
                }
            }   
        }
        if (partitioning_->ownPartitionContainsBottomBoundary() && !cellInfo.faceRight.isBoundaryFace()) {
            if (j == 2) {
                // bottom has to be a ghost cell. we have to set u there since paraview output writer needs that ghost u value to interpolate u on the boundary
                if (cellInfo.faceBottom.dirichletU.has_value()) {
                    double faceBCValue = cellInfo.faceBottom.dirichletU.value();
                    double u_i_j = discretization_->u(i,j);
                    discretization_->u(i,j-1) = 2.0 * faceBCValue - u_i_j;
                } else if (cellInfo.faceBottom.neumannU.has_value()) {
                    // set ghost value based on neumann condition
                    double neumannValue = cellInfo.faceBottom.neumannU.value();
                    double u_i_j = discretization_->u(i,j);
                    discretization_->u(i,j-1) = u_i_j + neumannValue * discretization_->dy();
                }
            }   
        }
        if (partitioning_->ownPartitionContainsRightBoundary()) {
            if (i == discretization_->nCells()[0] + 1) {
                // right has to be a ghost cell. we have to set v there since paraview output writer needs that ghost v value to interpolate v on the boundary
                if (cellInfo.faceRight.dirichletV.has_value()) {
                    double faceBCValue = cellInfo.faceRight.dirichletV.value();
                    double v_i_j = discretization_->v(i,j);
                    discretization_->v(i+1,j) = 2.0 * faceBCValue - v_i_j;
                } else if (cellInfo.faceRight.neumannV.has_value()) {
                    // set ghost value based on neumann condition
                    double neumannValue = cellInfo.faceRight.neumannV.value();
                    double v_i_j = discretization_->v(i,j);
                    discretization_->v(i+1,j) = v_i_j + neumannValue * discretization_->dx();
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