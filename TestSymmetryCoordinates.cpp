#include <cmath>

#include "SymmetryCoordinates.hpp"
#include "catch.hpp"

TEST_CASE("Cartesian Coordinate Assignment", "[cartcoord][short]") {
    SymmetryCoordinates::CartesianPoint2<double> p = {10.0, 242};
    REQUIRE(p.x == 10.0);
    REQUIRE(p.y == 242);
}

TEST_CASE("Symmetry Coordinate Assignment", "[symcoord][short]") {
    SymmetryCoordinates::SymmetryPoint2<double> sp = {10.0, 242};
    REQUIRE(sp.r == 10.0);
    REQUIRE(sp.t == 242);
}

TEST_CASE("Symmetry Coordinates From Cartesian", "[symcoord][short]") {
    SymmetryCoordinates::CartesianPoint2<double> cp0 = {0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2<double> cp1 = {1.0/3.0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2<double> cp2 = {-1.0/3.0, 1.0/std::sqrt(3)};
    SymmetryCoordinates::CartesianPoint2<double> cp3 = {2.0/3.0, 0};

    auto sp0 = SymmetryCoordinates::fromCartesian(cp0);
    auto sp1 = SymmetryCoordinates::fromCartesian(cp1);
    auto sp2 = SymmetryCoordinates::fromCartesian(cp2);
    auto sp3 = SymmetryCoordinates::fromCartesian(cp3);

    REQUIRE(sp0.r == Approx(1));
    REQUIRE(sp0.t == Approx(0));
    REQUIRE(sp1.r == Approx(1));
    REQUIRE(sp1.t == Approx(M_PI / 6));
    REQUIRE(sp2.r == Approx(1));
    REQUIRE(sp2.t == Approx(M_PI / 6));
    REQUIRE(sp3.r == Approx(1));
    REQUIRE(sp3.t == Approx(M_PI / 6));
}

TEST_CASE("Cartesian Coordinate Operators", "[cartcoord][short]") {
    SymmetryCoordinates::CartesianPoint2<double> p = {10.0, 242};
    SymmetryCoordinates::CartesianPoint2<double> p2 = {10.0, 242};
    SymmetryCoordinates::CartesianPoint2<double> q = {103.0, 12};

    REQUIRE(p == p2);
    REQUIRE(p != q);

    auto sum = p+q;
    REQUIRE(sum.x == Approx(113.0));
    REQUIRE(sum.y == Approx(254.0));

    auto diff = p-q;
    REQUIRE(diff.x == Approx(-93.0));
    REQUIRE(diff.y == Approx(230.0));
}
