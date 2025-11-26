#include "pressureSolver/RedBlackGaussSeidel.hpp"

RedBlackGaussSeidel::RedBlackGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning) {}

void RedBlackGaussSeidel::solve() {
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();
    const double lek = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps_2 = epsilon_ * epsilon_;

    int iter = 0;

    computeResidualNorm();
    
    while (iter < maximumNumberOfIterations_ && residualNorm_ > eps_2) {
        iter++;
        
        // Red-Black ordering: (i+j) % 2 determines color
        // Red cells: (i+j) % 2 == 0
        int optionalShift = 0;
        if ((partitioning_->nodeOffset()[0] + partitioning_->nodeOffset()[1]) % 2 == 1) {
            optionalShift = 1;
        }
            

        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
            const int i_start = discretization_->pIBegin() + ((j + optionalShift) % 2);
            for (int i = i_start; i <= discretization_->pIEnd(); i+=2) {
                    const double ersterTerm = (discretization_->p(i+1,j) + discretization_->p(i-1,j)) / dx2;
                    const double zweiterTerm = (discretization_->p(i,j+1) + discretization_->p(i,j-1)) / dy2;
                    discretization_->p(i,j) = lek * (ersterTerm + zweiterTerm - discretization_->rhs(i,j));
                }
            }
        
        // Communicate red cell values to neighbors
        communicateAndSetBoundaryValues();
        
        // Black cells: (i+j) % 2 == 1
        for (int j = discretization_->pJBegin(); j <= discretization_->pJEnd(); j++) {
            const int i_start = discretization_->pIBegin() + ((j + optionalShift + 1) % 2);
            for (int i = i_start; i <= discretization_->pIEnd(); i+=2) {
                    const double ersterTerm = (discretization_->p(i+1,j) + discretization_->p(i-1,j)) / dx2;
                    const double zweiterTerm = (discretization_->p(i,j+1) + discretization_->p(i,j-1)) / dy2;
                    discretization_->p(i,j) = lek * (ersterTerm + zweiterTerm - discretization_->rhs(i,j));
                }
            }
        
        // Communicate black cell values and set boundary values
        communicateAndSetBoundaryValues();
        
        // Compute residual norm after full iteration
        computeResidualNorm();
    }
    
    this->numberOfIterations_ = iter;
}