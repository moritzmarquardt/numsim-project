#include "parallelComputation.hpp"

void ParallelComputation::initialize(int argc, char *argv[]) {
    //TODO: implement parallel initialization
    settings_.loadFromFile(argv[1]);
    settings_.printSettings();

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
        // TODO: implement ParallelSOR
        std::cerr << "Error: Parallel SOR not yet implemented." << std::endl;
        std::exit(EXIT_FAILURE);
        // pressureSolver_ = std::make_unique<SOR>(discretization_, settings_.epsilon, settings_.maximumNumberOfIterations, settings_.omega);
    } else {
        std::cerr << "Error: Unknown pressure solver: " << settings_.pressureSolver << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // create output writers
    outputWriterParaview_ = std::make_unique<OutputWriterParaviewParallel>(discretization_, *partitioning_);
    outputWriterText_ = std::make_unique<OutputWriterTextParallel>(discretization_, *partitioning_);
}

void ParallelComputation::runSimulation() {
    //TODO: implement parallel simulation run
}

void ParallelComputation::computeTimeStepWidth() {
    // compute time step width based on CFL condition
    double dx = discretization_->dx();
    double dy = discretization_->dy();

    double dx2 = dx * dx;
    double dy2 = dy * dy;  

    double dt_diff_cond = 0.5 * settings_.re * (dx2 * dy2) / (dx2 + dy2);
    double dt_conv_cond_local = std::min(dx / discretization_->u().computeMaxAbs(), dy / discretization_->v().computeMaxAbs());

    double dt_conv_cond_global;
    MPI_Allreduce(&dt_conv_cond_local, &dt_conv_cond_global, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

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
    std::vector<double> sendBufferTopU(uIEnd - uIBegin + 1, 0.0);
    std::vector<double> sendBufferTopV(vIEnd - vIBegin + 1, 0.0);
    std::vector<double> sendBufferBottomU(uIEnd - uIBegin + 1, 0.0);
    std::vector<double> sendBufferBottomV(vIEnd - vIBegin + 1, 0.0);
    std::vector<double> sendBufferLeftU(uJEnd - uJBegin + 1, 0.0);
    std::vector<double> sendBufferLeftV(vJEnd - vJBegin + 1, 0.0);
    std::vector<double> sendBufferRightU(uJEnd - uJBegin + 1, 0.0);
    std::vector<double> sendBufferRightV(vJEnd - vJBegin + 1, 0.0);

    MPI_Request requestsTop, requestsBottom, requestsLeft, requestsRight;

    if (partitioning_->ownPartitionContainsTopBoundary()) {
        for (int i = uIBegin; i <= uIEnd; i++) { 
            discretization_->u(i,uJEnd + 1) = 2 * settings_.dirichletBcTop[0] - discretization_->u(i,uJEnd); 
        }
    } else {
        for (int i = uIBegin; i <= uIEnd; i++) {
            sendBufferTopU[i - uIBegin] = discretization_->u(i,uJEnd);
        }
        MPI_Isend(sendBufferTopU.data(), sendBufferTopU.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsTop);

        for (int i = vIBegin; i <= vIEnd; i++) {
            sendBufferTopV[i - vIBegin] = discretization_->v(i,vJEnd - 1);
        }
        MPI_Isend(sendBufferTopV.data(), sendBufferTopV.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 1, MPI_COMM_WORLD, &requestsTop);

        MPI_Irecv(sendBufferTopU.data(), sendBufferTopU.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 0, MPI_COMM_WORLD, &requestsTop);
        MPI_Irecv(sendBufferTopV.data(), sendBufferTopV.size(), MPI_DOUBLE, partitioning_->topNeighbourRankNo(), 1, MPI_COMM_WORLD, &requestsTop);
    }
    
    if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = uIBegin; i <= uIEnd; i++) { 
            discretization_->u(i,uJBegin - 1) = 2 * settings_.dirichletBcBottom[0] - discretization_->u(i,uJBegin); 
        }
    } else {
        for (int i = uIBegin; i <= uIEnd; i++) {
            sendBufferBottomU[i - uIBegin] = discretization_->u(i,uJBegin);
        }
        MPI_Isend(sendBufferBottomU.data(), sendBufferBottomU.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 2, MPI_COMM_WORLD, &requestsBottom);

        for (int i = vIBegin; i <= vIEnd; i++) {
            sendBufferBottomV[i - vIBegin] = discretization_->v(i,vJBegin);
        }
        MPI_Isend(sendBufferBottomV.data(), sendBufferBottomV.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 3, MPI_COMM_WORLD, &requestsBottom);

        MPI_Irecv(sendBufferBottomU.data(), sendBufferBottomU.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 2, MPI_COMM_WORLD, &requestsBottom);
        MPI_Irecv(sendBufferBottomV.data(), sendBufferBottomV.size(), MPI_DOUBLE, partitioning_->bottomNeighbourRankNo(), 3, MPI_COMM_WORLD, &requestsBottom);
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