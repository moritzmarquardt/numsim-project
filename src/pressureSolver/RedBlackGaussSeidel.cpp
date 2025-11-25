#include "pressureSolver/RedBlackGaussSeidel.hpp"

RedBlackGaussSeidel::RedBlackGaussSeidel(std::shared_ptr<Discretization> discretization, double epsilon, int maximumNumberOfIterations, std::shared_ptr<Partitioning> partitioning) :
    ParallelPressureSolver(discretization, epsilon, maximumNumberOfIterations, partitioning) {}

void RedBlackGaussSeidel::solve() {
    //TODO: implement Red-Black Gauss-Seidel method for parallel pressure solving
}