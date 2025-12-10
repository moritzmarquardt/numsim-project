#include "computation/parallelComputation.hpp"
#include <iostream>

#include <cstdlib>
#include "settings.hpp"
#include <mpi.h>
#include "storage/array2d.hpp"

// Helper function to create a vector of evenly spaced values
std::vector<double> linspace(double start, double end, int num) {
  std::vector<double> result(num);
  for (int i = 0; i < num; ++i) {
    result[i] = start + i * (end - start) / (num - 1);
  }
  return result;
}

int main(int argc, char *argv[])
{
  // emasure time for parallel execution
  MPI_Init(&argc, &argv);
  int rank;
  int num = 101;
  std::vector<double> Re = linspace(500, 1500, num);
  std::vector<double> vel = linspace(0.5, 1.5, num);
    
  ParallelComputation computation;
  computation.initialize(argc, argv);

  for (int i = 0; i < num; ++i) {
    computation.setRe(Re[i]);
    computation.setTopUBoundary(vel[i]);
    computation.setSimNumber(i);

    rank = computation.getRankNo();

    double startTime = MPI_Wtime();
    computation.runSimulation();
    double endTime = MPI_Wtime();
    if (rank == 0) {
      std::cout << "Run " << i+1 << "/" << num << " with Re = " << Re[i] << " and top boundary velocity = " << vel[i] << std::endl;
      std::cout << "Total simulation time: " << endTime - startTime << " seconds." << std::endl;
      std::cout << "----------------------------------------" << std::endl;
    }
  }

  MPI_Finalize();

  return EXIT_SUCCESS;
}