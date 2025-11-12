#pragma once

#include "array2d.hpp"
#include <array>
#include <iostream>

/** A field variable is the discretization of a scalar function f(x) with x in the computational domain. 
 * More specifically, a scalar value is stored at discrete nodes/points. 
 * The nodes are arranged in an equidistant mesh with specified mesh width. 
 */

class FieldVariable : public Array2D {
public:
    /**
     * Construct a field variable.

    Parameters:
        size	    The number of entries in x and y direction.
        origin	    Cartesian coordinates of the point with (i,j) = (0,0), this is different from (0,0) for the u,v and p field variables.
        meshWidth	The length of a single cell in x and y direction. 
     */
    FieldVariable(std::array<int,2> size, std::array<double,2> origin, std::array<double,2> meshWidth);

    FieldVariable& operator=(const FieldVariable& other);

    /**
     * get the value at the Cartesian coordinate (x,y). The value is linearly interpolated between stored points. 
     * We use bilinear interpolation to get the value at (x,y) based on the four surrounding points.
     */
    double interpolateAt(double x, double y) const;

    /**
     * Compute maximum absolute value for the field variable to use it in the time step condition.
     */
    double computeMaxAbs() const;

    /**
     * Print the field variable as 2D array to the standard output pretty formatted for debugging the code.
     */
    void printAsArray() const;

private:
    std::array<double,2> origin_; //Cartesian coordinates of the point with (i,j) = (0,0), this is different from (0,0) for the u,v and p field variables. 
    std::array<double,2> meshWidth_; //The length of a single cell in x and y direction. 
};
  
