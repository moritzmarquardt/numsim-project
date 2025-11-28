#include "computation/parallelComputation.hpp"
#include <iostream>

#include <cstdlib>
#include "settings.hpp"
#include <mpi.h>
#include "storage/array2d.hpp"

int main(int argc, char *argv[])
{

  if (argc < 2)
  {
    std::cerr << "Usage: " << (argc > 0 ? argv[0] : "numsim") << " <settings-file>" << std::endl;
    return EXIT_FAILURE;
  }
  else
  {
    // check if the settings file exists and can be opened
    std::ifstream settingsFile(argv[1]);
    if (!settingsFile.is_open())
    {
      std::cerr << "Error: Could not open settings file: " << argv[1] << std::endl;
      return EXIT_FAILURE;
    }
  }

  // emasure time for parallel execution
  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  ParallelComputation computation;
  computation.initialize(argc, argv);

  double startTime = MPI_Wtime();
  computation.runSimulation();
  double endTime = MPI_Wtime();

  MPI_Finalize();
  if (rank == 0) {
    std::cout << "Total simulation time: " << endTime - startTime << " seconds." << std::endl;
  }

  return EXIT_SUCCESS;
}