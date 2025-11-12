#include "FieldVariable.hpp"
#include <cassert>
#include <cmath>
#include <iostream>

// constructor with inherited Array2D part
FieldVariable::FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth) :
    Array2D(size), origin_(origin), meshWidth_(meshWidth) {}

/**
 * Implements the assignment operator for FieldVariable
 */
FieldVariable& FieldVariable::operator=(const FieldVariable& other) {
    if (this != &other) {
        Array2D::operator=(other); // Copy base class part
        origin_[0] = other.origin_[0];
        origin_[1] = other.origin_[1];
        meshWidth_[0] = other.meshWidth_[0];
        meshWidth_[1] = other.meshWidth_[1];
    }
    return *this;
}

double FieldVariable::interpolateAt(double x, double y) const {
    // Check if x and y are located in the domain
    assert(0.0 <= x && x < size_[0]*meshWidth_[0]);
    assert(0.0 <= y && y < size_[1]*meshWidth_[1]);

    // Find the indicies i and j of the cell which contains the left bottom corner of the square in which (x,y) is located.
    // origin is the cartesian coordinates of the field values in the (0,0) cell.
    // in order to get the index of the bottom left corner, we have to "revert" the shift of the origin
    const int i = static_cast<int>((x - origin_[0]) / meshWidth_[0]); 
    const int j = static_cast<int>((y - origin_[1]) / meshWidth_[1]);

    const double x_left = origin_[0] + i * meshWidth_[0];
    const double y_bottom = origin_[1] + j * meshWidth_[1];

    const double value_left_bottom = (*this)(i,j);
    const double value_right_bottom = (*this)(i + 1,j);
    const double value_left_top = (*this)(i, j + 1);
    const double value_right_top = (*this)(i + 1, j + 1);

    const double interpolate_bottom = value_left_bottom + (value_right_bottom - value_left_bottom) * (x - x_left) / meshWidth_[0];
    const double interpolate_top = value_left_top + (value_right_top - value_left_top) * (x - x_left) / meshWidth_[0];
    
    return interpolate_bottom + (interpolate_top - interpolate_bottom) * (y - y_bottom) / meshWidth_[1];
}

double FieldVariable::computeMaxAbs() const {
    double max_value = 0.0;

    // Go through all points
    for (int i = 0; i < size_[0]; i++) {
        for (int j = 0; j < size_[1]; j++) {
            double abs_val = fabs((*this)(i,j));
            if (abs_val > max_value) {
                max_value = abs_val;
            }
        }     
    }

    return max_value;
}

void FieldVariable::printAsArray() const {
    // pretty print the field variable as 2D array
    for (int j = size_[1] - 1; j >= 0; j--) {
        for (int i = 0; i < size_[0]; i++) {
            std::cout << (*this)(i, j) << " ";
        }
        std::cout << std::endl;
    }
}