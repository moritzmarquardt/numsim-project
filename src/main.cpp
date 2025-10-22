#include "output_writer/write_paraview_output.h"


#include <iostream>

#include <cstdlib>
#include "settings.h"


int main(int argc, char *argv[]){

  for (int i = 0; i < 5; i++)

  {

    writeParaviewOutput(i);

  }

  return EXIT_SUCCESS;

}