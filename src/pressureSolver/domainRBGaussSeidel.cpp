#include "pressureSolver/domainRBGaussSeidel.hpp"

DomainRBGaussSeidel::DomainRBGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning, std::shared_ptr<Domain> domain) :
    DomainPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning, domain) {}

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
                // if (partitioning_->ownPartitionContainsTopBoundary() && j == discretization_->pJEnd()) {
                //     // Dirichlet BC at top boundary
                //     discretization_->p(i, j + 1) = p_i_jp1; // assuming ghost cell already set to BC value
                // }
            }
            if (redCellInfo.bottomIsBoundaryFace()) {
                p_i_jm1 = p_i_j;
                // if (partitioning_->ownPartitionContainsBottomBoundary() && j == discretization_->pJBegin()) {
                //     // Dirichlet BC at bottom boundary
                //     discretization_->p(i, j - 1) = p_i_jm1; // assuming ghost cell already set to BC value
                // }
            }
            if (redCellInfo.leftIsBoundaryFace()) {
                p_im1_j = p_i_j;
                // if (partitioning_->ownPartitionContainsLeftBoundary() && i == discretization_->pIBegin()) {
                //     // Dirichlet BC at left boundary
                //     discretization_->p(i - 1, j) = p_im1_j; // assuming ghost cell already set to BC value
                // }
            }
            if (redCellInfo.rightIsBoundaryFace()) {
                p_ip1_j = p_i_j;
                // if (partitioning_->ownPartitionContainsRightBoundary() && i == discretization_->pIEnd()) {
                //     // Dirichlet BC at right boundary
                //     discretization_->p(i + 1, j) = p_ip1_j; // assuming ghost cell already set to BC value
                // }
            }

            discretization_->p(i,j) = lek * ((p_ip1_j + p_im1_j)/dx2 + (p_i_jp1 + p_i_jm1) / dy2 - discretization_->rhs(i,j));

            // // Todo: on partition boundary set Neumann zero BC after udating pressure.

            // if (redCellInfo.topIsBoundaryFace() && partitioning_->ownPartitionContainsTopBoundary() && j == discretization_->pJEnd()) {
            //     // Dirichlet BC at top boundary
            //     discretization_->p(i, j + 1) = discretization_->p(i, j); // assuming ghost cell already set to BC value
            // }
            // if (redCellInfo.bottomIsBoundaryFace() && partitioning_->ownPartitionContainsBottomBoundary() && j == discretization_->pJBegin()) {
            //     // Dirichlet BC at bottom boundary
            //     discretization_->p(i, j - 1) = discretization_->p(i, j); // assuming ghost cell already set to BC value
            // }
            // if (redCellInfo.leftIsBoundaryFace() && partitioning_->ownPartitionContainsLeftBoundary() && i == discretization_->pIBegin()) {
            //     // Dirichlet BC at left boundary
            //     discretization_->p(i - 1, j) = discretization_->p(i, j); // assuming ghost cell already set to BC value
            // }
            // if (redCellInfo.rightIsBoundaryFace() && partitioning_->ownPartitionContainsRightBoundary() && i == discretization_->pIEnd()) {
            //     // Dirichlet BC at right boundary
            //     discretization_->p(i + 1, j) = discretization_->p(i, j); // assuming ghost cell already set to BC value
            // }
        }

       
        // print debug info
        if (iter == 1 && partitioning_->ownRankNo() == 3) {
            std::cout << "After red update in iteration " << iter << ", p field:" << std::endl;
            for (int j = discretization_->pJEnd()+1; j >= discretization_->pJBegin()-1; j--) {
                for (int i = discretization_->pIBegin()-1; i <= discretization_->pIEnd()+1; i++) {
                    std::cout << discretization_->p(i, j) << " ";
                }
                std::cout << std::endl;
            }
        }
        

         // debug shit
        const int pIBegin = discretization_->pIBegin();
        const int pIEnd = discretization_->pIEnd();
        const int pJBegin = discretization_->pJBegin();
        const int pJEnd = discretization_->pJEnd();
        if (partitioning_->ownPartitionContainsTopBoundary()) {
            for (int i = pIBegin; i <= pIEnd; i++) {
                discretization_->p(i, pJEnd + 1) = discretization_->p(i, pJEnd);
            }
        }
        if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            discretization_->p(i, pJBegin - 1) = discretization_->p(i, pJBegin);
        }
        }
        if (partitioning_->ownPartitionContainsLeftBoundary()) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                discretization_->p(pIBegin - 1, j) = discretization_->p(pIBegin, j);
            }
        }
        if (partitioning_->ownPartitionContainsRightBoundary()) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                discretization_->p(pIEnd + 1, j) = discretization_->p(pIEnd, j);
            }
        }
        
        // Communicate red cell values to neighbors
        communicateGhostValues();
        if (iter == 1 && partitioning_->ownRankNo() == 3) {
            std::cout << "After communicating red cells in iteration " << iter << ", p field:" << std::endl;
            for (int j = discretization_->pJEnd()+1; j >= discretization_->pJBegin()-1; j--) {
                for (int i = discretization_->pIBegin()-1; i <= discretization_->pIEnd()+1; i++) {
                    std::cout << discretization_->p(i, j) << " ";
                }
                std::cout << std::endl;
            }
        }


        std::vector<CellInfo> blackCellsInfo = domain_->getBlackListFluid();
        n = blackCellsInfo.size();
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
                // if (partitioning_->ownPartitionContainsTopBoundary() && j == discretization_->pJEnd()) {
                //     // Dirichlet BC at top boundary
                //     discretization_->p(i, j + 1) = p_i_jp1; // assuming ghost cell already set to BC value
                // }
            }
            if (blackCellInfo.bottomIsBoundaryFace()) {
                p_i_jm1 = p_i_j;
                // if (partitioning_->ownPartitionContainsBottomBoundary() && j == discretization_->pJBegin()) {
                //     // Dirichlet BC at bottom boundary
                //     discretization_->p(i, j - 1) = p_i_jm1; // assuming ghost cell already set to BC value
                // }
            }
            if (blackCellInfo.leftIsBoundaryFace()) {
                p_im1_j = p_i_j;
                // if (partitioning_->ownPartitionContainsLeftBoundary() && i == discretization_->pIBegin()) {
                //     // Dirichlet BC at left boundary
                //     discretization_->p(i - 1, j) = p_im1_j; // assuming ghost cell already set to BC value
                // }
            }
            if (blackCellInfo.rightIsBoundaryFace()) {
                p_ip1_j = p_i_j;
                // if (partitioning_->ownPartitionContainsRightBoundary() && i == discretization_->pIEnd()) {
                //     // Dirichlet BC at right boundary
                //     discretization_->p(i + 1, j) = p_ip1_j; // assuming ghost cell already set to BC value
                // }
            }

            discretization_->p(i,j) = lek * ((p_ip1_j + p_im1_j)/dx2 + (p_i_jp1 + p_i_jm1) / dy2 - discretization_->rhs(i,j));

            if (blackCellInfo.topIsBoundaryFace() && partitioning_->ownPartitionContainsTopBoundary() && j == discretization_->pJEnd()) {
                // Dirichlet BC at top boundary
                discretization_->p(i, j + 1) = discretization_->p(i, j); // assuming ghost cell already set to BC value
            }
            if (blackCellInfo.bottomIsBoundaryFace() && partitioning_->ownPartitionContainsBottomBoundary() && j == discretization_->pJBegin()) {
                // Dirichlet BC at bottom boundary
                discretization_->p(i, j - 1) = discretization_->p(i, j); // assuming ghost cell already set to BC value
            }
            if (blackCellInfo.leftIsBoundaryFace() && partitioning_->ownPartitionContainsLeftBoundary() && i == discretization_->pIBegin()) {
                // Dirichlet BC at left boundary
                discretization_->p(i - 1, j) = discretization_->p(i, j); // assuming ghost cell already set to BC value
            }
            if (blackCellInfo.rightIsBoundaryFace() && partitioning_->ownPartitionContainsRightBoundary() && i == discretization_->pIEnd()) {
                // Dirichlet BC at right boundary
                discretization_->p(i + 1, j) = discretization_->p(i, j); // assuming ghost cell already set to BC value
            }
        }

        // set boundary values before communicating black cells
        if (partitioning_->ownPartitionContainsTopBoundary()) {
            for (int i = pIBegin; i <= pIEnd; i++) {
                discretization_->p(i, pJEnd + 1) = discretization_->p(i, pJEnd);
            }
        }
        if (partitioning_->ownPartitionContainsBottomBoundary()) {
        for (int i = pIBegin; i <= pIEnd; i++) {
            discretization_->p(i, pJBegin - 1) = discretization_->p(i, pJBegin);
        }
        }
        if (partitioning_->ownPartitionContainsLeftBoundary()) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                discretization_->p(pIBegin - 1, j) = discretization_->p(pIBegin, j);
            }
        }
        if (partitioning_->ownPartitionContainsRightBoundary()) {
            for (int j = pJBegin; j <= pJEnd; j++) {
                discretization_->p(pIEnd + 1, j) = discretization_->p(pIEnd, j);
            }
        }

        
        // Communicate black cell values and set boundary values
        communicateGhostValues();
        
        // Compute residual norm after full iteration
        computeResidualNorm();
    }
    
    this->numberOfIterations_ = iter;
}