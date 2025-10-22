#include "output_writer/write_paraview_output.hpp"


#include <iostream>

#include <cstdlib>
#include "settings.hpp"
#include "storage/array2d.hpp"


int main(int argc, char *argv[]){
  
  Array2D myArray = Array2D({2, 2});

  myArray(0,0) = 7.0;
  myArray(1,0) = 2.0;
  std::cout << "Array values:" << myArray(0,0) << "marquis penis" << myArray(1,1)<< std::endl;

  std::array<int,2> size = myArray.size();

  return EXIT_SUCCESS;

}