#include "output_writer/write_paraview_output.h"


#include <iostream>

#include <cstdlib>
#include "settings.h"


int main(int argc, char *argv[]){

  Settings settings;
  settings.loadFromFile("parameters.txt");
  settings.printSettings();

  return EXIT_SUCCESS;

}