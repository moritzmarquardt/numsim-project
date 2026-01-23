#include "pressureSolver/domainRBGaussSeidel.hpp"

DomainRBGaussSeidel::DomainRBGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    DomainPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning) {}

void DomainRBGaussSeidel::solve() {
    const double dx2 = discretization_->dx() * discretization_->dx();
    const double dy2 = discretization_->dy() * discretization_->dy();
    const double lek = (dx2 * dy2) / (2.0 * (dx2 + dy2));
    const double eps_2 = epsilon_ * epsilon_;

    int iter = 0;

    computeResidualNorm();
    
    while (iter < maximumNumberOfIterations_ && residualNorm_ > eps_2) {
        iter++;

        std::vector<CellInfo> redCellsInfo = domain_->getRedListFluid();
        int n = redCellsInfo.size();
        double dx = discretization_->dx();
        for (int idx = 0; idx < n; idx++) {
            CellInfo redCellInfo = redCellsInfo[idx];
            int i = redCellInfo.cellIndexPartition[0];
            int j = redCellInfo.cellIndexPartition[1];

            double p_i_j = discretization_->p(i,j);
            double p_ip1_j = discretization_->p(i+1,j);
            double p_im1_j = discretization_->p(i-1,j);
            double p_i_jp1 = discretization_->p(i,j+1);
            double p_i_jm1 = discretization_->p(i,j-1);

            if (redCellInfo.topIsBoundaryFace()) {
                p_i_jp1 = p_i_j;
            }
            if (redCellInfo.bottomIsBoundaryFace()) {
                p_i_jm1 = p_i_j;
            }
            if (redCellInfo.leftIsBoundaryFaceC()) {
                p_im1_j = p_i_j;
            }
            if (redCellInfo.rightIsBoundaryFace()) {
                p_ip1_j = p_i_j;
            }

            discretization_->p(i,j) = lek * ((p_ip1_j + p_im1_j)/dx2 + (p_i_jp1 + p_i_jm1) / dy2 - discretization_->rhs(i,j));
        }
        
        // Communicate red cell values to neighbors
        communicateAndSetBoundaryValues();


        std::vector<CellInfo> blackCellsInfo = domain_->getBlackListFluid();
        int n = blackCellsInfo.size();
        double dx = discretization_->dx();
        for (int idx = 0; idx < n; idx++) {
            CellInfo blackCellInfo = blackCellsInfo[idx];
            int i = blackCellInfo.cellIndexPartition[0];
            int j = blackCellInfo.cellIndexPartition[1];

            double p_i_j = discretization_->p(i,j);
            double p_ip1_j = discretization_->p(i+1,j);
            double p_im1_j = discretization_->p(i-1,j);
            double p_i_jp1 = discretization_->p(i,j+1);
            double p_i_jm1 = discretization_->p(i,j-1);

            if (blackCellInfo.topIsBoundaryFace()) {
                p_i_jp1 = p_i_j;
            }
            if (blackCellInfo.bottomIsBoundaryFace()) {
                p_i_jm1 = p_i_j;
            }
            if (blackCellInfo.leftIsBoundaryFaceC()) {
                p_im1_j = p_i_j;
            }
            if (blackCellInfo.rightIsBoundaryFace()) {
                p_ip1_j = p_i_j;
            }

            discretization_->p(i,j) = lek * ((p_ip1_j + p_im1_j)/dx2 + (p_i_jp1 + p_i_jm1) / dy2 - discretization_->rhs(i,j));
        }
        
        // Communicate black cell values and set boundary values
        communicateAndSetBoundaryValues();
        
        // Compute residual norm after full iteration
        computeResidualNorm();
    }
    
    this->numberOfIterations_ = iter;
}