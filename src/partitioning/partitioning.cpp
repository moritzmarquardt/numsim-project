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
    const int periods[2] = {0, 0}; // non-periodic
    MPI_Cart_create(MPI_COMM_WORLD, dimensions, nSubdomains_.data(), periods, 1, &cartComm_);

    // get own coordinates in Cartesian communicator
    MPI_Comm_rank(cartComm_, &ownRankNo_);                        // set ownRankNo_
    MPI_Cart_coords(cartComm_, ownRankNo_, dimensions, ownCoords_.data()); // set ownCoords_

    // calculate how many cells are assigned to each subdomain
    for (int dim = 0; dim < dimensions; dim++)
    {
        std::array<int,2> dist = calcDist1D(nCellsGlobal_[dim], nSubdomains_[dim], ownCoords_[dim]);
        nCellsLocal_[dim] = dist[0];
        nodeOffset_[dim] = dist[1];
    }

    // calculate neighbouring rank numbers
    MPI_Cart_shift(cartComm_, 0, 1, &leftNeighbourRankNo_, &rightNeighbourRankNo_);  // shift in dimension 0 (x-direction)
    MPI_Cart_shift(cartComm_, 1, 1, &bottomNeighbourRankNo_, &topNeighbourRankNo_);  // shift in dimension 1 (y-direction)

    std::cout << "\n=== Rank " << ownRankNo_
            << " / " << nRanks_-1 << " ===\n";
    std::cout << "Local cells: [" << nCellsLocal_[0] 
                << ", " << nCellsLocal_[1] << "]\n";
    std::cout << "Node offset: [" << nodeOffset_[0] 
                << ", " << nodeOffset_[1] << "]\n";
    std::cout << "Boundaries: "
                << (ownPartitionContainsLeftBoundary() ? "LEFT " : "")
                << (ownPartitionContainsRightBoundary() ? "RIGHT " : "")
                << (ownPartitionContainsBottomBoundary() ? "BOTTOM " : "")
                << (ownPartitionContainsTopBoundary() ? "TOP " : "")
                << "\n";
    std::cout << "Rank " << ownRankNo_ << " has coordinates (" << ownCoords_[0] << "," << ownCoords_[1] << ")\n";
    std::cout << "Rank has neighbours: left " << leftNeighbourRankNo_ << ", right " << rightNeighbourRankNo_
              << ", bottom " << bottomNeighbourRankNo_ << ", top " << topNeighbourRankNo_ << "\n";


}

std::array<int,2> Partitioning::nCellsLocal() const
{
    return nCellsLocal_;
}

std::array<int,2> Partitioning::nCellsGlobal() const
{
    return nCellsGlobal_;
}

int Partitioning::ownRankNo() const
{
    return ownRankNo_;
}

int Partitioning::nRanks() const
{
    return nRanks_;
}

bool Partitioning::ownPartitionContainsBottomBoundary() const
{
    return (ownCoords_[1] == 0);
}

bool Partitioning::ownPartitionContainsTopBoundary() const
{
    return (ownCoords_[1] == nSubdomains_[1] - 1);
}

bool Partitioning::ownPartitionContainsLeftBoundary() const
{
    return (ownCoords_[0] == 0);
}

bool Partitioning::ownPartitionContainsRightBoundary() const
{
    return (ownCoords_[0] == nSubdomains_[0] - 1);
}

int Partitioning::leftNeighbourRankNo() const
{
    return leftNeighbourRankNo_;
}

int Partitioning::rightNeighbourRankNo() const
{
    return rightNeighbourRankNo_;
}

int Partitioning::topNeighbourRankNo() const
{
    return topNeighbourRankNo_;
}

int Partitioning::bottomNeighbourRankNo() const
{
    return bottomNeighbourRankNo_;
}

std::array<int,2> Partitioning::nodeOffset() const
{
    return nodeOffset_;
}