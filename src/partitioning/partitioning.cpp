#include "partitioning.hpp"
#include <mpi.h>
#include <iostream>

std::array<int,2> Partitioning::calcDist1D(int nCellsGlobal, int nSubdomains, int ownCoord)
{
    std::array<int,2> result{0,0}; // {nCellsLocal, nodeOffset}

    // cells that can be evenly distributed
    const int nCellsDist = nCellsGlobal / nSubdomains; // implicit int casting

    // cells that cannot be evenly distributed
    const int nCellsRemaining = nCellsGlobal % nSubdomains;

    if (ownCoord < nCellsRemaining)
    {
        result[0] = nCellsDist + 1;
        result[1] = ownCoord * (nCellsDist + 1);
    }
    else
    {
        result[0] = nCellsDist;
        result[1] = nCellsRemaining * (nCellsDist + 1) + (ownCoord - nCellsRemaining) * nCellsDist;
    }

    return result;
}

void Partitioning::initialize(std::array<int, 2> nCellsGlobal)
{
    nCellsGlobal_ = nCellsGlobal;
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks_); // set nRanks_

    // let MPI compute a good division of ranks in x and y direction
    nSubdomains_ = {0, 0};
    const int dimensions = 2;
    MPI_Dims_create(nRanks_, dimensions, nSubdomains_.data()); // set nSubdomains_

    // create Cartesian communicator based on the computed number of subdomains
    MPI_Comm cartComm;
    const int periods[2] = {0, 0}; // non-periodic
    MPI_Cart_create(MPI_COMM_WORLD, dimensions, nSubdomains_.data(), periods, 0, &cartComm);

    // get own coordinates in Cartesian communicator
    MPI_Comm_rank(cartComm, &ownRankNo_);                        // set ownRankNo_
    MPI_Cart_coords(cartComm, ownRankNo_, dimensions, ownCoords_.data()); // set ownCoords_
    std::cout << "Rank " << ownRankNo_ << " has coordinates (" << ownCoords_[0] << "," << ownCoords_[1] << ")\n";

    // calculate how many cells are assigned to each subdomain
    for (int dim = 0; dim < dimensions; dim++)
    {
        std::array<int,2> dist = calcDist1D(nCellsGlobal_[dim], nSubdomains_[dim], ownCoords_[dim]);
        nCellsLocal_[dim] = dist[0];
        nodeOffset_[dim] = dist[1];
    }

}