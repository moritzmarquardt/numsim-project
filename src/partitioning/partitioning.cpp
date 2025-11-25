#include "partitioning.hpp"
#include <mpi.h>

void Partitioning::initialize(std::array<int,2> nCellsGlobal)
{
    nSubdomains_ = {0,0};
    nCellsGlobal_ = nCellsGlobal;

    MPI_Dims_create(nRanks_, 2, nSubdomains_.data());
    // create Cartesian communicator
    MPI_Comm cartComm;
    const int periods[2] = {0, 0}; // non-periodic
    MPI_Cart_create(MPI_COMM_WORLD, 2, nSubdomains_.data(), periods, 0, &cartComm);

    // get own coordinates in Cartesian communicator
    MPI_Comm_rank(cartComm, &ownRankNo_);
    MPI_Comm_size(cartComm, &nRanks_);
    MPI_Cart_coords(cartComm, ownRankNo_, 2, ownCoords_.data());

    const int nCellsRemainingX = nCellsGlobal_[0] % nSubdomains_[0]; // left over cells when division is not perfect
    const int nCellsRemainingY = nCellsGlobal_[1] % nSubdomains_[1];

    if (ownCoords_[0] < nCellsRemainingX)
    {
        // TODO: check how ranks are ordered and the indexing of ownCoords_
        nCellsLocal_[0] = nCellsGlobal_[0] / nSubdomains_[0] + 1;
        nodeOffset_[0] = ownCoords_[0] * nCellsLocal_[0];
    }
    else
    {
        nCellsLocal_[0] = nCellsGlobal_[0] / nSubdomains_[0];
        nodeOffset_[0] = nCellsRemainingX * (nCellsGlobal_[0] / nSubdomains_[0] + 1)
                        + (ownCoords_[0] - nCellsRemainingX) * nCellsLocal_[0];
    }


    
}


/*
std::array<int,2> nCellsLocal_;    //< number of cells in own partition
    std::array<int,2> nodeOffset_;     //< offset of nodes in own partition
    */