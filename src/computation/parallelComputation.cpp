#include "parallelComputation.hpp"
#include <cmath>

void ParallelComputation::initialize(int argc, char *argv[]) {
    //TODO: implement parallel initialization
    settings_.loadFromFile(argv[1]);
    // settings_.printSettings();

    partitioning_ = std::make_shared<Partitioning>();
    partitioning_->initialize(settings_.nCells);
    std::array<int,2> nCellsLocal = partitioning_->nCellsLocal();

    // calculate mesh width
    meshWidth_[0] = settings_.physicalSize[0] / settings_.nCells[0];
    meshWidth_[1] = settings_.physicalSize[1] / settings_.nCells[1];

    // create discretization
    if (settings_.useDonorCell) {
        discretization_ = std::make_shared<DonorCell>(nCellsLocal, meshWidth_, settings_.alpha, partitioning_);
    } else {
        discretization_ = std::make_shared<CentralDifferences>(nCellsLocal, meshWidth_,partitioning_);
    }

    // std::cout << "Created discretization with mesh width dx: " << meshWidth_[0] << ", dy: " << meshWidth_[1] << std::endl;
    // std::cout << "Number of cells in x direction: " << nCellsLocal[0] << ", y direction: " << nCellsLocal[1] << std::endl;

    // create pressure solver
    if (settings_.pressureSolver == "GaussSeidel") {
        pressureSolver_ = std::make_unique<RedBlackGaussSeidel>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, partitioning_);
    } else if (settings_.pressureSolver == "SOR") {
        // TODO: calculate optimal omega (not in settings hardcoded)
        pressureSolver_ = std::make_unique<RedBlackSOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega, partitioning_);
    } else {
        std::cerr << "Error: Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // create output writers
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, *partitioning_);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, *partitioning_);
}

void ParallelComputation::runSimulation() {
    applyInitialBoundaryValues();

    double currentTime = 0.0;
    int iterationCount = 0;
    const double time_eps = 1e-8;

    while (currentTime < settings_.endTime - time_eps) {
        applyBoundaryValues();

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
        // if (iterationCount % 10 == 0 || currentTime >= settings_.endTime && partitioning_->ownRankNo() == 0) {
        //     int percent = static_cast<int>((currentTime / settings_.endTime) * 100);
            
        //     // Create progress bar
        //     const int barWidth = 40;
        //     int pos = barWidth * currentTime / settings_.endTime;
        //     std::string progressBar = "[";
        //     for (int i = 0; i < barWidth; ++i) {
        //         if (i < pos) progressBar += "=";
        //         else if (i == pos) progressBar += ">";
        //         else progressBar += " ";
        //     }
        //     progressBar += "]";
            
        //     std::cout << "\rProgress: " << progressBar << " " << percent << "% | Time: " << currentTime 
        //             << "/" << settings_.endTime << " | Iter: " << iterationCount << std::flush;
        // }
        // std::cout << "Advanced to time: " << currentTime << std::endl;
        // std::cout << "Completed iteration: " << iterationCount << std::endl;

        outputWriterParaview_->writeFile(currentTime);
        // std::cout << "Wrote Paraview output." << std::endl;
        // outputWriterText_->writeFile(currentTime);
        // std::cout << "Wrote text output." << std::endl;

        // std::cout << "Iteration: " << iterationCount << ", Time: " << currentTime << ", dt: " << dt_ << std::endl;
        
    }
    // std::cout << std::endl << "Simulation completed." << std::endl;
}

void ParallelComputation::computeTimeStepWidth() {
    // compute time step width based on CFL condition
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    double dx2 = dx * dx;
    double dy2 = dy * dy;  

    double dt_diff_cond = 0.5 * settings_.re * (dx2 * dy2) / (dx2 + dy2);
    double dt_conv_cond_local = std::min(dx / discretization_->u().computeMaxAbs(), dy / discretization_->v().computeMaxAbs());

    MPI_Request time_request;
    double dt_conv_cond_global = 0.0;
    MPI_Iallreduce(&dt_conv_cond_local, &dt_conv_cond_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD, &time_request);
    MPI_Wait(&time_request, MPI_STATUS_IGNORE);

    double dt_prelim = settings_.tau * std::min(dt_diff_cond, dt_conv_cond_global);

    if (dt_prelim < settings_.maximumDt) {
        dt_ = dt_prelim;
    } else {
        dt_ = settings_.maximumDt;
        std::cout << "Warning: Time step width limited by maximumDt!" << std::endl;
    }
}

void ParallelComputation::applyBoundaryValues() {
    const int uIBegin = discretization_->uIBegin();
    const int uIEnd = discretization_->uIEnd();
    const int uJBegin = discretization_->uJBegin();
    const int uJEnd = discretization_->uJEnd();

    const int vIBegin = discretization_->vIBegin();
    const int vIEnd = discretization_->vIEnd();
    const int vJBegin = discretization_->vJBegin();
    const int vJEnd = discretization_->vJEnd();

    // buffers
    std::vector<double> sendBufferTop(uIEnd - uIBegin + vIEnd - vIBegin + 2, 0.0);
    std::vector<double> sendBufferBottom(uIEnd - uIBegin + vIEnd - vIBegin + 2, 0.0);
    std::vector<double> sendBufferLeft(uJEnd - uJBegin + vJEnd - vJBegin + 2, 0.0);
    std::vector<double> sendBufferRight(uJEnd - uJBegin + vJEnd - vJBegin + 2, 0.0);

    MPI_Request requestsTop, requestsBottom, requestsLeft, requestsRight;
    const int TAG_U = 0;
    const int TAG_V = 1;

    if (partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = uIBegin; i <= uIEnd; i++) { 
            discretization_->u(i,uJEnd + 1) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i,uJEnd); 
        }
    } else {
        for (int i = uIBegin; i <= uIEnd; i++) {
            sendBufferTop[i - uIBegin] = discretization_->u(i,uJEnd);
        }

        for (int i = vIBegin; i <= vIEnd; i++) {
            sendBufferTop[i - vIBegin + uIEnd + 1 - uIBegin] = discretization_->v(i,vJEnd - 1);
        }
        MPI_Isend(sendBufferTop.data(), sendBufferTop.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsTop);

        MPI_Irecv(sendBufferTop.data(), sendBufferTop.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsTop);
    }
    
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = uIBegin; i <= uIEnd; i++) { 
            discretization_->u(i,uJBegin - 1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i,uJBegin); 
        }
    } else {
        for (int i = uIBegin; i <= uIEnd; i++) {
            sendBufferBottom[i - uIBegin] = discretization_->u(i,uJBegin);
        }

        for (int i = vIBegin; i <= vIEnd; i++) {
            sendBufferBottom[i - vIBegin + uIEnd + 1 - uIBegin] = discretization_->v(i,vJBegin + 1); // +1 because we have two layers of gjost cells at the bottom (just like at the left)
        }
        MPI_Isend(sendBufferBottom.data(), sendBufferBottom.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsBottom);

        MPI_Irecv(sendBufferBottom.data(), sendBufferBottom.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), TAG_U, MPI_COMM_WORLD, &requestsBottom);
    }

    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = vJBegin - 1; j <= vJEnd + 1; j++) { 
            discretization_->v(vIBegin - 1,j) = 2 * settings_.dirichletBcLeft[1] - discretization_->v(vIBegin,j); 
        }
    } else {
        for (int j = uJBegin; j <= uJEnd; j++) {
            sendBufferLeft[j - uJBegin] = discretization_->u(uIBegin + 1,j); // +1 because we have two layers of ghost cells at the left (just like at the bottom)
        }

        for (int j = vJBegin; j <= vJEnd; j++) {
            sendBufferLeft[j - vJBegin + uJEnd + 1 - uJBegin] = discretization_->v(vIBegin,j);
        }
        MPI_Isend(sendBufferLeft.data(), sendBufferLeft.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsLeft);
        
        MPI_Irecv(sendBufferLeft.data(), sendBufferLeft.size(), MPI_DOUBLE, partitioning_->leftNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsLeft);
    }

    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = vJBegin - 1; j <= vJEnd + 1; j++) { 
            discretization_->v(vIEnd + 1,j) = 2 * settings_.dirichletBcRight[1] - discretization_->v(vIEnd,j); 
        }
    } else {
        for (int j = uJBegin; j <= uJEnd; j++) {
            sendBufferRight[j - uJBegin] = discretization_->u(uIEnd - 1,j);
        }

        for (int j = vJBegin; j <= vJEnd; j++) {
            sendBufferRight[j - vJBegin + uJEnd + 1 - uJBegin] = discretization_->v(vIEnd,j);
        }
        MPI_Isend(sendBufferRight.data(), sendBufferRight.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsRight);
        
        MPI_Irecv(sendBufferRight.data(), sendBufferRight.size(), MPI_DOUBLE, partitioning_->rightNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsRight);
    }

    // wait for all communications to finish
    if (!partitioning_->ownPartitionContainsTopBoundary()) {
        MPI_Wait(&requestsTop, MPI_STATUS_IGNORE);
        for (int i = uIBegin; i <= uIEnd; i++) {
            discretization_->u(i,uJEnd + 1) = sendBufferTop[i - uIBegin];
        }
        for (int i = vIBegin; i <= vIEnd; i++) {
            discretization_->v(i,vJEnd + 1) = sendBufferTop[i - vIBegin + uIEnd + 1 - uIBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsBottomBoundary()) {
        MPI_Wait(&requestsBottom, MPI_STATUS_IGNORE);
        for (int i = uIBegin; i <= uIEnd; i++) {
            discretization_->u(i,uJBegin - 1) = sendBufferBottom[i - uIBegin];
        }
        for (int i = vIBegin; i <= vIEnd; i++) {
            discretization_->v(i,vJBegin - 1) = sendBufferBottom[i - vIBegin + uIEnd + 1 - uIBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsLeftBoundary()) {
        MPI_Wait(&requestsLeft, MPI_STATUS_IGNORE);
        for (int j = uJBegin; j <= uJEnd; j++) {
            discretization_->u(uIBegin - 1,j) = sendBufferLeft[j - uJBegin];
        }
        for (int j = vJBegin; j <= vJEnd; j++) {
            discretization_->v(vIBegin - 1,j) = sendBufferLeft[j - vJBegin + uJEnd + 1 - uJBegin];
        }
    }
    if (!partitioning_->ownPartitionContainsRightBoundary()) {
        MPI_Wait(&requestsRight, MPI_STATUS_IGNORE);
        for (int j = uJBegin; j <= uJEnd; j++) {
            discretization_->u(uIEnd + 1,j) = sendBufferRight[j - uJBegin];
        }
        for (int j = vJBegin; j <= vJEnd; j++) {
            discretization_->v(vIEnd + 1,j) = sendBufferRight[j - vJBegin + uJEnd + 1 - uJBegin];
        }
    }
}

void ParallelComputation::applyInitialBoundaryValues() {
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++) {
            const int j_bottom_bc = discretization_->vJBegin() - 1;
            discretization_->v(i, j_bottom_bc) = settings_.dirichletBcBottom[1]; // bottom
            discretization_->g(i, j_bottom_bc) = settings_.dirichletBcBottom[1]; // bottom
        }
    }
    if (partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = discretization_->vIBegin(); i <= discretization_->vIEnd(); i++) {
            const int j_top_bc = discretization_->vJEnd() + 1;
            discretization_->v(i, j_top_bc) = settings_.dirichletBcTop[1]; // top
            discretization_->g(i, j_top_bc) = settings_.dirichletBcTop[1]; // top
        }
    }
    if (partitioning_->ownPartitionContainsLeftBoundary()) {
        for (int j = discretization_->uJBegin()-1; j <= discretization_->uJEnd()+1; j++) {
            const int i_left_bc = discretization_->uIBegin() - 1;
            discretization_->u(i_left_bc, j) = settings_.dirichletBcLeft[0]; // left
            discretization_->f(i_left_bc, j) = settings_.dirichletBcLeft[0]; // left
        }
    }
    if (partitioning_->ownPartitionContainsRightBoundary()) {
        for (int j = discretization_->uJBegin()-1; j <= discretization_->uJEnd()+1; j++) {
            const int i_right_bc = discretization_->uIEnd() + 1;
            discretization_->u(i_right_bc, j) = settings_.dirichletBcRight[0]; // right
            discretization_->f(i_right_bc, j) = settings_.dirichletBcRight[0]; // right
        }
    }
}