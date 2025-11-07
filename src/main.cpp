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

  Array2D myArray = Array2D({2, 2});

  myArray(0, 0) = 7.0;
  myArray(1, 0) = 2.0;
  std::cout << "Array values (selected): " << myArray(0, 0) << " marquise " << myArray(1, 1) << std::endl;

  std::array<int, 2> size = myArray.size();

  // iterate over all elements and print them
  for (int j = 0; j < size[1]; ++j)
  {
    for (int i = 0; i < size[0]; ++i)
    {
      std::cout << "myArray(" << i << "," << j << ") = " << myArray(i, j) << std::endl;
      myArray(i, j) = i + j;
      std::cout << "myArray(" << i << "," << j << ") = " << myArray(i, j) << std::endl;
    }
  }

  return EXIT_SUCCESS;
}