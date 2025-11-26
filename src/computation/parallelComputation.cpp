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
    //TODO: implement parallel time step width computation
}

void ParallelComputation::applyBoundaryValues() {
    //TODO: implement parallel boundary value application
}

void ParallelComputation::applyInitialBoundaryValues() {
    //TODO: implement parallel initial boundary value application
}