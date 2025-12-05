#include "computation/parallelComputation.hpp"
#include <iostream>

#include <cstdlib>
#include "settings.hpp"
#include <mpi.h>
#include "storage/array2d.hpp"

int main(int argc, char *argv[])
{
  // emasure time for parallel execution
  MPI_Init(&argc, &argv);
  int rank;

  ParallelComputation computation;
  computation.initialize(argc, argv);
  rank = computation.getRankNo();

  double startTime = MPI_Wtime();
  computation.runSimulation();
  double endTime = MPI_Wtime();

  MPI_Finalize();
  if (rank == 0) {
    std::cout << "Total simulation time: " << endTime - startTime << " seconds." << std::endl;
  }

  return EXIT_SUCCESS;
}