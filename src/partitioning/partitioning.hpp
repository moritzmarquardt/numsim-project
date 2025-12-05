#pragma once

#include <array>
#include <mpi.h>

/**
 * @class Partitioning
 * @brief Manages domain decomposition for parallel computation using MPI.
 * 
 * Attributes:
 * - nCellsLocal_: Number of cells in the local subdomain.
 * - nCellsGlobal_: Total number of cells in the global domain.
 * - nodeOffset_: Offset of nodes in the local subdomain relative to the global domain.
 * - ownRankNo_: The MPI rank number of the current process.
 * - nRanks_: Total number of MPI ranks.
 * - nSubdomains_: Number of subdomains in each spatial direction.
 */
class Partitioning
{
public:

  //! compute partitioning, set internal variables
  void initialize(std::array<int,2> nCellsGlobal);

  // calculate distribution of nodes along one axis
  // distributes the given number of global cells to the given number of subdomains
  // important if the number of cells is not divisible by the number of subdomains
  // return array with {nCellsLocal, nodeOffset}
  std::array<int,2> calcDist1D(int nCellsGlobal, int nSubdomains, int ownCoord);

  //! get the local number of cells in the own subdomain
  std::array<int,2> nCellsLocal() const;

  //! get the global number of cells in the whole computational domain
  //! used in OutputWriterParaviewParallel
  std::array<int,2> nCellsGlobal() const;

  //! get the own MPI rank no
  //! used in OutputWriterParaviewParallel and OutputWriterTextParallel
  int ownRankNo() const;

  //! number of MPI ranks
  int nRanks() const;

  //! if the own partition has part of the bottom boundary of the whole domain
  bool ownPartitionContainsBottomBoundary() const;

  //! if the own partition has part of the top boundary of the whole domain
  //! used in OutputWriterParaviewParallel
  bool ownPartitionContainsTopBoundary() const;

  //! if the own partition has part of the left boundary of the whole domain
  bool ownPartitionContainsLeftBoundary() const;

  //! if the own partition has part of the right boundary of the whole domain
  //! used in OutputWriterParaviewParallel
  bool ownPartitionContainsRightBoundary() const;

  //! get the rank no of the left neighbouring rank
  int leftNeighbourRankNo() const;

  //! get the rank no of the right neighbouring rank
  int rightNeighbourRankNo() const;

  //! get the rank no of the top neighbouring rank
  int topNeighbourRankNo() const;

  //! get the rank no of the bottom neighbouring rank
  int bottomNeighbourRankNo() const;

  //! get the offset values for counting local nodes in x and y direction. 
  //! (i_local,j_local) + nodeOffset = (i_global,j_global)
  //! used in OutputWriterParaviewParallel
  std::array<int,2> nodeOffset() const;

  //! get the global number of cells in the whole computational domain
  int getNCellsGlobal() const;

  //! get the Cartesian communicator
  MPI_Comm getCartComm() const;

  private:
    std::array<int,2> nCellsLocal_;    //< number of cells in own partition
    std::array<int,2> nCellsGlobal_;   //< global number of cells
    std::array<int,2> nodeOffset_;     //< offset of nodes in own partition
    int ownRankNo_;                    //< own MPI rank no
    std::array<int,2> ownCoords_;               //< own coordinates in the Cartesian communicator
    int nRanks_;                       //< number of MPI ranks
    std::array<int,2> nSubdomains_;   //< number of subdomains in x and y direction
    MPI_Comm cartComm_;              //< Cartesian communicator
    int leftNeighbourRankNo_;      //< rank no of left neighbour
    int rightNeighbourRankNo_;     //< rank no of right neighbour
    int topNeighbourRankNo_;       //< rank no of top neighbour
    int bottomNeighbourRankNo_;    //< rank no of bottom neighbour
};
