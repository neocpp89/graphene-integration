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

TEST_CASE("Interpolation Evaluation (w2)", "[interpolations][short]") {
    const int Nr = 100;
    const int Nt = 50;
    const double dr = 1.0 / Nr;
    const double dt = 2 * M_PI / Nt;

    for(int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            const double r = i * dr;
            const double t = j * dt;
            CHECK(w2_naive(r, t) == Approx(w2_horner(r, t)));
        }
    }
}

TEST_CASE("Interpolation Evaluation (w3)", "[interpolations][short]") {
    const int Nr = 100;
    const int Nt = 50;
    const double dr = 1.0 / Nr;
    const double dt = 2 * M_PI / Nt;

    for(int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            const double r = i * dr;
            const double t = j * dt;
            CHECK(w3_naive(r, t) == Approx(w3_horner(r, t)));
        }
    }
}

TEST_CASE("Interpolation Evaluation (w4)", "[interpolations][short]") {
    const int Nr = 100;
    const int Nt = 50;
    const double dr = 1.0 / Nr;
    const double dt = 2 * M_PI / Nt;

    for(int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            const double r = i * dr;
            const double t = j * dt;
            CHECK(w4_naive(r, t) == Approx(w4_horner(r, t)));
        }
    }
}

TEST_CASE("Interpolation Evaluation (w5)", "[interpolations][short]") {
    const int Nr = 100;
    const int Nt = 50;
    const double dr = 1.0 / Nr;
    const double dt = 2 * M_PI / Nt;

    for(int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            const double r = i * dr;
            const double t = j * dt;
            CHECK(w5_naive(r, t) == Approx(w5_horner(r, t)));
        }
    }
}

TEST_CASE("Interpolation Evaluation (w6)", "[interpolations][short]") {
    const int Nr = 100;
    const int Nt = 50;
    const double dr = 1.0 / Nr;
    const double dt = 2 * M_PI / Nt;

    for(int i = 0; i < Nr; i++) {
        for (int j = 0; j < Nt; j++) {
            const double r = i * dr;
            const double t = j * dt;
            CHECK(w6_naive(r, t) == Approx(w6_horner(r, t)));
        }
    }
}
