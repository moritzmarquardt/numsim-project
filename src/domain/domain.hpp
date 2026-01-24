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

        /**
         * Read domain file and set up the domain including all the lists of cells
         */
        void readDomainFile(const std::string& filename);

        // getter for boundaryInfoListAll_
        std::vector<CellInfo> getInfoListAll() const {
            return *cellListAll_;
        }

        /**
         * Get all cells with info that are fluid cells (no obstacle cells)
         */
        std::vector<CellInfo> getInfoListFluid() const {
            return *cellListFluid_;
        }
        std::vector<CellInfo> getRedListFluid() const {
            return *redListFluid_;
        }
        std::vector<CellInfo> getBlackListFluid() const {
            return *blackListFluid_;
        }

    private:
        const Settings* settings_;
        std::shared_ptr<Partitioning> partitioning_;
        std::unique_ptr<Array2D> obstacleMask_; // has size of local partition
        std::unique_ptr<std::vector<CellInfo>> cellListAll_; // list of all cells in the local partition with their cell info
        std::unique_ptr<std::vector<CellInfo>> cellListFluid_; // list of all fluid cells in the local partition with their cell info
        std::unique_ptr<std::vector<CellInfo>> redListFluid_; // list of all red fluid cells in the local partition with their cell info
        std::unique_ptr<std::vector<CellInfo>> blackListFluid_; // list of all black fluid cells in the local partition with their cell info
};