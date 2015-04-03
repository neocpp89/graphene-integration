#include <iostream>
#include <cmath>
#include "SymmetryCoordinates.hpp"

int main()
{
    SymmetryCoordinates::CartesianPoint2 cp0 = {0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2 cp1 = {1.0/3.0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2 cp2 = {-1.0/3.0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2 cp3 = {2.0/3.0, 0};

    SymmetryCoordinates::SymmetryPoint2 sp0 = SymmetryCoordinates::fromCartesian(cp0);
    SymmetryCoordinates::SymmetryPoint2 sp1 = SymmetryCoordinates::fromCartesian(cp1);
    SymmetryCoordinates::SymmetryPoint2 sp2 = SymmetryCoordinates::fromCartesian(cp2);
    SymmetryCoordinates::SymmetryPoint2 sp3 = SymmetryCoordinates::fromCartesian(cp3);

    std::cout << sp0.r << ", " << sp0.t << '\n';
    std::cout << sp1.r << ", " << sp1.t << '\n';
    std::cout << sp2.r << ", " << sp2.t << '\n';
    std::cout << sp3.r << ", " << sp3.t << '\n';

    return 0;
}
