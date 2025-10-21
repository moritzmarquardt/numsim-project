#include "settings.h"

#include <fstream> // for file operations
#include <iostream> // for cout

/**
 * Parse a text file with settings, each line contains "<parameterName> = <value>"
 * 
 * @param filename Name of the file to load settings from
 */
void Settings::loadFromFile(std::string filename){
  // open file
  std::ifstream file(filename.c_str(), std::ios::in);

  // check if file is open
  if (!file.is_open()){
    std::cout << "Could not open parameter file \"" << filename << "\"." << std::endl;
    return;
  }

  // loop over lines of file
  for (int lineNo = 0;; lineNo++){
    // read line
    std::string line;
    getline(file, line);

    // at the end of the file break for loop
    if (file.eof())
      break;

    // print line 
    std::cout << "line " << lineNo << ": " << line << std::endl;
  }
}

/**
 * Output all settings to console
 */
void Settings::printSettings(){
  std::cout << "Settings: " << std::endl
            << "  physicalSize: " << physicalSize[0] << " x " << physicalSize[1] << ", nCells: " << nCells[0] << " x " << nCells[1] << std::endl
            << "  endTime: " << endTime << " s, re: " << re << ", g: (" << g[0] << "," << g[1] << "), tau: " << tau << ", maximum dt: " << maximumDt << std::endl
            << "  dirichletBC: bottom: (" << dirichletBcBottom[0] << "," << dirichletBcBottom[1] << ")"
            << ", top: (" << dirichletBcTop[0] << "," << dirichletBcTop[1] << ")"
            << ", left: (" << dirichletBcLeft[0] << "," << dirichletBcLeft[1] << ")"
            << ", right: (" << dirichletBcRight[0] << "," << dirichletBcRight[1] << ")" << std::endl
            << "  useDonorCell: " << std::boolalpha << useDonorCell << ", alpha: " << alpha << std::endl
            << "  pressureSolver: " << pressureSolver << ", omega: " << omega << ", epsilon: " << epsilon << ", maximumNumberOfIterations: " << maximumNumberOfIterations << std::endl;
}