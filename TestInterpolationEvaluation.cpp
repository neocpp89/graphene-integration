#include <cmath>

#include "Interpolations.hpp"
#include "catch.hpp"

TEST_CASE("Interpolation Evaluation (w1)", "[interpolations][short]") {
    const int Nr = 100;
    const int Nt = 50;
    const double dr = 1.0 / Nr;
    const double dt = 2 * M_PI / Nt;

    for(int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            const double r = i * dr;
            const double t = j * dt;
            CHECK(w1_naive(r, t) == Approx(w1_horner(r, t)));
        }
    }
}

