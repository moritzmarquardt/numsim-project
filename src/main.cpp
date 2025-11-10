#include "output_writer/write_paraview_output.hpp"
#include "computation/computation.hpp"
#include <iostream>

#include <cstdlib>
#include "settings.hpp"
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

  

  Computation computation;
  std::cout << "Init..." << std::endl;
  computation.initialize(argc, argv);
  std::cout << "Initialization done. Running simulation..." << std::endl;
  computation.runSimulation();

  return EXIT_SUCCESS;
}