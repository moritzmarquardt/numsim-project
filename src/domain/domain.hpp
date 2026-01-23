#pragma once

#include <array>
#include <memory>
#include <mpi.h>
#include "storage/array2d.hpp"
#include "discretization/staggeredGrid.hpp"
#include "partitioning/partitioning.hpp"
#include "settings.hpp"
#include "discretization/discretization.hpp"
#include "domain/boundaryInfo.hpp"
#include "domain/cellInfo.hpp"

class Domain {
    public:
        Domain(const Settings* settings, std::shared_ptr<Partitioning> partitioning);

        void readDomainFile(const std::string& filename);

        /**
         * Get the fluid mask for the local partition.
         */
        std::unique_ptr<Array2D> fluidMaskPartition();

        // getter for boundaryInfoListAll_
        std::vector<CellInfo> getBoundaryInfoListAll() {
            return *boundaryInfoListAll_;
        }

    private:
        const Settings* settings_;
        std::shared_ptr<Partitioning> partitioning_;
        std::unique_ptr<Array2D> fluidMask_; 
        std::unique_ptr<Array2D> obstacleMask_;
        std::unique_ptr<std::vector<CellInfo>> boundaryInfoListAll_; //holds all cells with their respective boundary condition info. All cells means including the additional "ghost" row and column left and below the partition which store the info for the outer face of the partition.
        // std::unique_ptr<std::vector<std::pair<int, int>>> uIndexes_; // list of (i,j) indexes where u has to be calculated (equivalent to ui begin und end und uj begin und end functions in StaggeredGrid) so basically fluid cells that are inner cells.
        std::unique_ptr<std::vector<BoundaryInfo>> uList; // list of boundary infos for all boundaries in the domain
        // std::unique_ptr<std::vector<std::pair<int, int>>> vIndexes_; // list of (i,j) indexes where v has to be calculated
        std::unique_ptr<std::vector<BoundaryInfo>> vList; // list of boundary infos for all boundaries in the domain
        // std::unique_ptr<std::vector<std::pair<int, int>>> pIndexes_; // list of (i,j) indexes where p has to be calculated (this is equivalent to the list of fluid cells since p is in the middle of the cell and never is on the same spot as a boundary conditions bc they are only defined on faces /edges)
        std::unique_ptr<std::vector<BoundaryInfo>> pList; // list of boundary infos for all boundaries in the domain
};