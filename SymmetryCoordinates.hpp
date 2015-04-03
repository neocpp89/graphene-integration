#ifndef SYMMETRYCOORDINATES
#define SYMMETRYCOORDINATES
#include <cmath>

namespace SymmetryCoordinates {

struct CartesianPoint2 {
    double x;
    double y;

    friend bool operator==(const CartesianPoint2 &a, const CartesianPoint2 &b);
    friend bool operator<(const CartesianPoint2 &a, const CartesianPoint2 &b);
    friend CartesianPoint2 operator-(const CartesianPoint2 &a, const CartesianPoint2 &b);
    friend CartesianPoint2 operator+(const CartesianPoint2 &a, const CartesianPoint2 &b);

};

bool operator==(const CartesianPoint2 &a, const CartesianPoint2 &b) {
    return (a.x == b.x && a.y == b.y);
}

bool operator<(const CartesianPoint2 &a, const CartesianPoint2 &b) {
    if (a.x == b.x) {
        return (a.y < b.y);
    } else {
        return (a.x < b.x);
    }
}

CartesianPoint2 operator-(const CartesianPoint2 &a, const CartesianPoint2 &b)
{
    CartesianPoint2 cp = {a.x-b.x, a.y-b.y};
    return cp;
}

CartesianPoint2 operator+(const CartesianPoint2 &a, const CartesianPoint2 &b)
{
    CartesianPoint2 cp = {a.x+b.y, a.y+b.y};
    return cp;
}

struct SymmetryPoint2 {
    double r;
    double t;
};

constexpr double L = std::sqrt(3.0);

double calculateTheta(const CartesianPoint2 &cp);
double calculateRho(const CartesianPoint2 &cp, const double theta);
double calculateRho(const CartesianPoint2 &cp);

SymmetryPoint2 fromCartesian(const CartesianPoint2 &cp)
{
    SymmetryPoint2 sp;
    sp.t = calculateTheta(cp);
    sp.r = calculateRho(cp, sp.t);
    return sp;
}
    
double calculateTheta(const CartesianPoint2 &cp)
{
    const double magcp = std::sqrt(cp.x*cp.x + cp.y*cp.y);
    if (magcp == 0) {
        return 0;
    }

    const double theta = std::acos(std::cos(6.0*(std::acos(cp.x / magcp) + M_PI / 2.0))) / 6.0;
    return theta;
}

double calculateRho(const CartesianPoint2 &cp, const double theta)
{
    const double magcp = std::sqrt(cp.x*cp.x + cp.y*cp.y);
    const double rho = L * magcp * std::cos(theta);
    return rho;
}

double calculateRho(const CartesianPoint2 &cp)
{
    return calculateRho(cp, std::cos(calculateTheta(cp)));
}


};

#endif //SYMMETRYCOORDINATES
