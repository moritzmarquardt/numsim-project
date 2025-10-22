#include "output_writer/write_paraview_output.hpp"


#include <iostream>

#include <cstdlib>
#include "settings.hpp"


int main(int argc, char *argv[]){

  for (int i = 0; i < 5; i++)

  {

    writeParaviewOutput(i);

  }

  

  std::cout << "Program finished successfully." << std::endl;

  return EXIT_SUCCESS;

}