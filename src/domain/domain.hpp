#pragma once

#include <array>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>
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
        /** 
         * @brief Construct a new Domain object
         * 
         * @param settings Pointer to the Settings object
         * @param partitioning Shared pointer to the Partitioning object
         */
        Domain(const Settings* settings, std::shared_ptr<Partitioning> partitioning);

        /**
         * @brief Read domain file and set up the domain including all the lists of cells
         * @param filename Name of the domain file
         */
        void readDomainFile(const std::string& filename);

        /**
         * @brief Create CellInfo for a given cell based on the arrays storing the obstacle mask and boundary conditions
         */
        CellInfo createCellInfo(int iGlobal, int jGlobal, int iLocal, int jLocal, int nCellsXLocal, int nCellsYLocal);

        // getter for boundaryInfoListAll_
        std::vector<CellInfo> getInfoListAll() const {
            return *cellListAllLocal_;
        }

        /**
         * @brief Get all cells with info that are fluid cells (no obstacle cells)
         */
        std::vector<CellInfo> getInfoListFluid() const {
            return *cellListFluidLocal_;
        }
        /**
         * @brief Get all red fluid cells with info
         */
        std::vector<CellInfo> getRedListFluid() const {
            return *redListFluidLocal_;
        }
        /**
         * @brief Get all black fluid cells with info
         */
        std::vector<CellInfo> getBlackListFluid() const {
            return *blackListFluidLocal_;
        }
        /**
         * @brief Get all ghost cells with info (cells adjacent to left and bottom partition boundaries)
         * This includes only the cells that are outside of the domain.
         */
        std::vector<CellInfo> getGhostList() const {
            return *ghostListLocal_;
        }

        const std::unordered_map<double, BoundaryInfo>& rightFaceBCInfoMap() const {
            return rightFaceBCInfo_;
        }

        const std::unordered_map<double, BoundaryInfo>& topFaceBCInfoMap() const {
            return topFaceBCInfo_;
        }

        const std::unordered_map<double, char>& rightFaceMarkerMap() const {
            return rightFaceCodeToMarker_;
        }

        const std::unordered_map<double, char>& topFaceMarkerMap() const {
            return topFaceCodeToMarker_;
        }

        const Settings* settings_;
        std::shared_ptr<Partitioning> partitioning_;
        std::unique_ptr<Array2D> obstacleMaskGlobal_; // has size of global partition
        std::unique_ptr<Array2D> rightFacesBCGlobal_; // has size of global  + 1 in x direction
        std::unique_ptr<Array2D> topFacesBCGlobal_; // has size of global partition + 1 in y direction
        std::unique_ptr<std::vector<CellInfo>> cellListAllLocal_; // list of all cells in the local partition with their cell info
        std::unique_ptr<std::vector<CellInfo>> cellListFluidLocal_; // list of all fluid cells in the local partition with their cell info
        std::unique_ptr<std::vector<CellInfo>> redListFluidLocal_; // list of all red fluid cells in the local partition with their cell info
        std::unique_ptr<std::vector<CellInfo>> blackListFluidLocal_; // list of all black fluid cells in the local partition with their cell info
        std::unique_ptr<std::vector<CellInfo>> ghostListLocal_; // list with the cells left and bottom ghost cells in the local partition with their cell info

        std::unordered_map<double, BoundaryInfo> rightFaceBCInfo_;
        std::unordered_map<double, BoundaryInfo> topFaceBCInfo_;
        std::unordered_map<double, char> rightFaceCodeToMarker_;
        std::unordered_map<double, char> topFaceCodeToMarker_;
};