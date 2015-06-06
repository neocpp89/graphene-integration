#ifndef SYMMETRYCOORDINATES
#define SYMMETRYCOORDINATES
#include <iostream>
#include <cmath>

namespace SymmetryCoordinates {

template <typename Real = double>
struct CartesianPoint2 {
    Real x;
    Real y;

    bool operator==(const CartesianPoint2 &b)
    {
        return (x == b.x && y == b.y);
    }

    bool operator!=(const CartesianPoint2 &b)
    {
        return !(*this == b);
    }

    CartesianPoint2 &operator-=(const CartesianPoint2 &b)
    {
        auto &a = *this;
        a.x -= b.x;
        a.y -= b.y;
        return a;
    }

    CartesianPoint2 &operator+=(const CartesianPoint2 &b)
    {
        auto &a = *this;
        a.x += b.x;
        a.y += b.y;
        return a;
    }

    CartesianPoint2 operator-(const CartesianPoint2 &b) const
    {
        auto a = *this;
        return (a -= b);
    }

    CartesianPoint2 operator+(const CartesianPoint2 &b) const
    {
        auto a = *this;
        return (a += b);
    }

    bool operator<(const CartesianPoint2 &b) const
    {
        auto &a = *this;
        if (a.x == b.x) {
            return (a.y < b.y);
        } else {
            return (a.x < b.x);
        }
    }
};

template <typename Real>
std::ostream &operator<<(std::ostream &os, const CartesianPoint2<Real> &p)
{
    os << p.x << ", " << p.y;
    return os;
}

template <typename Real = double>
struct SymmetryPoint2 {
    Real r;
    Real t;
};

const double L = std::sqrt(3.0);

template <typename Real>
Real calculateTheta(const CartesianPoint2<Real> &cp);
template <typename Real>
Real calculateRho(const CartesianPoint2<Real> &cp, const Real theta);
template <typename Real>
Real calculateRho(const CartesianPoint2<Real> &cp);

template <typename Real>
SymmetryPoint2<Real> fromCartesian(const CartesianPoint2<Real> &cp)
{
    SymmetryPoint2<Real> sp;
    sp.t = calculateTheta<Real>(cp);
    sp.r = calculateRho<Real>(cp, sp.t);
    return sp;
}
    
template <typename Real>
Real calculateTheta(const CartesianPoint2<Real> &cp)
{
    const Real magcp = sqrt(cp.x*cp.x + cp.y*cp.y);
    if (magcp == 0) {
        return 0;
    }

    const Real theta= acos(cos(6.0*(acos(cp.x / magcp) + M_PI / 2.0))) / 6.0;
    return theta;
}

template <typename Real>
Real calculateRho(const CartesianPoint2<Real> &cp, const Real theta)
{
    const Real magcp = sqrt(cp.x*cp.x + cp.y*cp.y);
    const Real rho = L * magcp * cos(theta);
    return rho;
}

template <typename Real>
Real calculateRho(const CartesianPoint2<Real> &cp)
{
    return calculateRho(cp, calculateTheta(cp));
}


};

#endif //SYMMETRYCOORDINATES
