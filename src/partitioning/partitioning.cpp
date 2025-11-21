#include "partitioning.hpp"
#include <mpi.h>

void Partitioning::initialize(std::array<int,2> nCellsGlobal)
{
    nSubdomains_ = {0,0};
    nCellsGlobal_ = nCellsGlobal;

    // get own rank no and number of ranks
    MPI_Comm_rank(MPI_COMM_WORLD, &ownRankNo_);
    MPI_Comm_size(MPI_COMM_WORLD, &nRanks_);
    MPI_Dims_create(nRanks_, 2, nSubdomains_.data());

    
}
