#include <cmath>

#include "SymmetryCoordinates.hpp"
#include "catch.hpp"

TEST_CASE("Cartesian Coordinate Assignment", "[cartcoord][short]") {
    SymmetryCoordinates::CartesianPoint2 p = {10.0, 242};
    REQUIRE(p.x == 10.0);
    REQUIRE(p.y == 242);
}

TEST_CASE("Symmetry Coordinate Assignment", "[symcoord][short]") {
    SymmetryCoordinates::SymmetryPoint2 sp = {10.0, 242};
    REQUIRE(sp.r == 10.0);
    REQUIRE(sp.t == 242);
}

TEST_CASE("Symmetry Coordinates From Cartesian", "[symcoord][short]") {
    SymmetryCoordinates::CartesianPoint2 cp0 = {0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2 cp1 = {1.0/3.0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2 cp2 = {-1.0/3.0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2 cp3 = {2.0/3.0, 0};

    SymmetryCoordinates::SymmetryPoint2 sp0 = SymmetryCoordinates::fromCartesian(cp0);
    SymmetryCoordinates::SymmetryPoint2 sp1 = SymmetryCoordinates::fromCartesian(cp1);
    SymmetryCoordinates::SymmetryPoint2 sp2 = SymmetryCoordinates::fromCartesian(cp2);
    SymmetryCoordinates::SymmetryPoint2 sp3 = SymmetryCoordinates::fromCartesian(cp3);

    REQUIRE(sp0.r == Approx(1));
    REQUIRE(sp0.t == Approx(0));
    REQUIRE(sp1.r == Approx(1));
    REQUIRE(sp1.t == Approx(M_PI / 6));
    REQUIRE(sp2.r == Approx(1));
    REQUIRE(sp2.t == Approx(M_PI / 6));
    REQUIRE(sp3.r == Approx(1));
    REQUIRE(sp3.t == Approx(M_PI / 6));
}
