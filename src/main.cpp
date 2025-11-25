#include "computation/computation.hpp"
#include "partitioning/partitioning.hpp"
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

  MPI_Init(&argc, &argv);

  Partitioning partitioning;
  partitioning.initialize({100, 100}); // Example global cell counts

  MPI_Barrier(MPI_COMM_WORLD);  // Synchronize output


  MPI_Finalize();

  return EXIT_SUCCESS;
}