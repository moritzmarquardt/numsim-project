#include "partitioning.hpp"
#include <mpi.h>
#include <algorithm>

void Partitioning::initialize(std::array<int,2> nCellsGlobal)
{
    nSubdomains_ = {0,0};
    nCellsGlobal_ = nCellsGlobal;

    // get own rank no and number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo_);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks_);
    MPI_Dims_create(nRanks_, 2, nSubdomains_.data());

    // Compute own coordinates in the 2D grid of subdomains.
    // MPI uses row-major ordering: rank = coords[1] * nSubdomains_[0] + coords[0]
    // So: coords[0] = rank % nSubdomains_[0]  (x-coordinate)
    //     coords[1] = rank / nSubdomains_[0]  (y-coordinate)
    ownCoords_[0] = ownRankNo_ % nSubdomains_[0];
    ownCoords_[1] = ownRankNo_ / nSubdomains_[0];

    // Calculate number of cells in local subdomain
    // Base cells per subdomain is nCellsGlobal / nSubdomains
    // Remainder cells are distributed to the first few subdomains
    int baseCellsX = nCellsGlobal_[0] / nSubdomains_[0];
    int remainderX = nCellsGlobal_[0] % nSubdomains_[0];
    int baseCellsY = nCellsGlobal_[1] / nSubdomains_[1];
    int remainderY = nCellsGlobal_[1] % nSubdomains_[1];

    // Subdomains with coords < remainder get one extra cell
    nCellsLocal_[0] = baseCellsX + (ownCoords_[0] < remainderX ? 1 : 0);
    nCellsLocal_[1] = baseCellsY + (ownCoords_[1] < remainderY ? 1 : 0);

    // Calculate node offset (starting position in global grid)
    // Subdomains before us have either baseCells+1 (if coord < remainder) or baseCells cells
    nodeOffset_[0] = ownCoords_[0] * baseCellsX + std::min(ownCoords_[0], remainderX);
    nodeOffset_[1] = ownCoords_[1] * baseCellsY + std::min(ownCoords_[1], remainderY);
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
    return ownCoords_[1] == 0;
}

bool Partitioning::ownPartitionContainsTopBoundary() const
{
    return ownCoords_[1] == nSubdomains_[1] - 1;
}

bool Partitioning::ownPartitionContainsLeftBoundary() const
{
    return ownCoords_[0] == 0;
}

bool Partitioning::ownPartitionContainsRightBoundary() const
{
    return ownCoords_[0] == nSubdomains_[0] - 1;
}

int Partitioning::leftNeighbourRankNo() const
{
    // Return -1 if there is no left neighbour (we are at left boundary)
    if (ownPartitionContainsLeftBoundary())
        return -1;
    
    // Left neighbour has coords[0] - 1 in x-direction, same y-coordinate
    // rank = coords[1] * nSubdomains_[0] + coords[0]
    return ownCoords_[1] * nSubdomains_[0] + (ownCoords_[0] - 1);
}

int Partitioning::rightNeighbourRankNo() const
{
    // Return -1 if there is no right neighbour (we are at right boundary)
    if (ownPartitionContainsRightBoundary())
        return -1;
    
    // Right neighbour has coords[0] + 1 in x-direction, same y-coordinate
    return ownCoords_[1] * nSubdomains_[0] + (ownCoords_[0] + 1);
}

int Partitioning::topNeighbourRankNo() const
{
    // Return -1 if there is no top neighbour (we are at top boundary)
    if (ownPartitionContainsTopBoundary())
        return -1;
    
    // Top neighbour has coords[1] + 1 in y-direction, same x-coordinate
    return (ownCoords_[1] + 1) * nSubdomains_[0] + ownCoords_[0];
}

int Partitioning::bottomNeighbourRankNo() const
{
    // Return -1 if there is no bottom neighbour (we are at bottom boundary)
    if (ownPartitionContainsBottomBoundary())
        return -1;
    
    // Bottom neighbour has coords[1] - 1 in y-direction, same x-coordinate
    return (ownCoords_[1] - 1) * nSubdomains_[0] + ownCoords_[0];
}

std::array<int,2> Partitioning::nodeOffset() const
{
    return nodeOffset_;
}
