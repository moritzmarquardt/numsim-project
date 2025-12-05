#include "array2d.hpp"

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