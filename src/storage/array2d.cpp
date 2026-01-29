#include "array2d.hpp"

#include <iostream>
#include <iomanip>
#include <stdexcept>

// constructor
Array2D::Array2D(std::array<int,2> size) :
  size_(size)
{
  // allocate data, initialize to 0
  data_.resize(size_[0]*size_[1], 0.0);
}

/**
 * Set all values to zero.
 */
void Array2D::setToZero() {
  std::fill(data_.begin(), data_.end(), 0.0);
}

/**
 * Get pointer to raw data.
 */
void* Array2D::data() {
  return data_.data();
}

/**
 * Assignment operator =
 */
Array2D& Array2D::operator=(const Array2D& other) {
    // Check if the sizes match
    if (size_ != other.size_) {
        throw std::runtime_error("Cannot assign Array2D objects of different sizes.");
    }

    // Copy the data
    data_ = other.data_;

    // Return a reference to this
    return *this;
}

void Array2D::printArray2D() {
    for (int j = this->size()[1] - 1; j >= 0; --j) {
        for (int i = 0; i < this->size()[0]; ++i) {
            std::cout << (*this)(i, j) << " ";
        }
        std::cout << std::endl;
    }
}

void Array2D::prettyPrintArray2D() {
    for (int j = this->size()[1] - 1; j >= 0; --j) {
        for (int i = 0; i < this->size()[0]; ++i) {
            std::cout << std::setw(8) << std::setprecision(4) << (*this)(i, j) << " ";
        }
        std::cout << std::endl;
    }
}