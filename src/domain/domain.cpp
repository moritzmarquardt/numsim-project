#include "domain.hpp"
#include <fstream>
#include <sstream>
#include <iostream>

Domain::Domain(const Settings* settings, std::shared_ptr<Partitioning> partitioning)
    : settings_(settings), partitioning_(partitioning) {
    // Constructor implementation is empty, nothing has to be done here. the interesting stuff happens in readDomainFile
}

void Domain::readDomainFile(const std::string& filename) {
    // Implement domain file reading logic here
    // This should read the domain configuration from the file
    std::ifstream file(filename);
    if (!file.is_open()) {
        std::cerr << "Error: Could not open domain file: " << filename << std::endl;
        return;
    }

    double physicalSizeX = settings_->physicalSize[0];
    double physicalSizeY = settings_->physicalSize[1];
    int nCellsX = settings_->nCells[0];
    int nCellsY = settings_->nCells[1];
    
    // first go through the whole file once to save grid info in arrays2d and marker info in a lsit of maps
    

    
    file.close();
}


