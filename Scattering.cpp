#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>

#include "fadiff.h"
#include "SymmetryCoordinates.hpp"
#include "PhysicalConstants.hpp"
#include "Interpolations.hpp"

const size_t POINTS_HINT = 4000000;

double delta(double arg, double sigma=10.0)
{
    return std::exp(arg*arg/(2*sigma*sigma));
}

template <typename Real>
Real PlanckDistribution(Real x)
{
    return 1.0 / (exp(x) - 1.0);
}

template <typename Real>
Real NSplitting(Real w_prime, Real w_doubleprime, Real temperature = 300)
{
    const Real factor = DIRAC_CONSTANT / (BOLTZMANN_CONSTANT * temperature);
    const Real x_prime = factor*w_prime;
    const Real x_doubleprime = factor*w_doubleprime;
    return (1 + PlanckDistribution(x_prime) + PlanckDistribution(x_doubleprime));
}

template <typename Real>
Real NCombining(Real w_prime, Real w_doubleprime, Real temperature = 300)
{
    const Real factor = DIRAC_CONSTANT / (BOLTZMANN_CONSTANT * temperature);
    const Real x_prime = factor*w_prime;
    const Real x_doubleprime = factor*w_doubleprime;
    return (PlanckDistribution(x_prime) - PlanckDistribution(x_doubleprime));
}

template <bool splitting, typename Real>
std::vector<SymmetryCoordinates::CartesianPoint2<Real>> CalculateQPrime(const SymmetryCoordinates::CartesianPoint2<Real> &q, const std::vector<SymmetryCoordinates::CartesianPoint2<Real>> &coord)
{
    std::vector<SymmetryCoordinates::CartesianPoint2<Real>> curve;
    curve.reserve(coord.size());

    if (splitting) {
        for (auto const &cp : coord) {
            auto diff = q - cp;
            auto sp = SymmetryCoordinates::fromCartesian(diff);
            if (sp.r > 1) { 
                curve.push_back(cp);
            }
        }
    } else {
        //combining process
        for (auto const &cp : coord) {
            auto diff = q + cp;
            auto sp = SymmetryCoordinates::fromCartesian(diff);
            if (sp.r > 1) { 
                curve.push_back(cp);
            }
        }
    }

    /* curve contains q' in cartesian coordinates. */
    return curve;
}

template <typename Real>
void domegadq(Real &xcomp, Real &ycomp, const SymmetryCoordinates::CartesianPoint2<Real> &q, const SymmetryCoordinates::SymmetryPoint2<Real> &qsym, fadbad::F<Real> (*fadfunc)(fadbad::F<Real>, fadbad::F<Real>))
{
    fadbad::F<double> sr,st;
    sr = qsym.r;
    st = qsym.t;
    sr.diff(0,2);
    st.diff(1,2);

    auto w = LATTICE_A*1e12*fadfunc(sr, st);
    const double dwdrg = w.d(0);
    const double dwdthetag = w.d(1);
    const double polar_r = sqrt(q.x*q.x + q.y*q.y);
    const double polar_t = atan2(q.y, q.x);
    const double thetag = qsym.t;
    const double modk = polar_r;
    const double dwdr = dwdrg * cos(thetag);
    const double dwdt = dwdrg * (-modk * sin(thetag)) + dwdthetag;
    const double dwdx = dwdr * cos(polar_t) - polar_r * dwdt * sin(polar_t);
    const double dwdy = dwdr * sin(polar_t) + polar_r * dwdt * cos(polar_t);
    xcomp = dwdx;
    ycomp = dwdy;
    return;
}

// q is the point to evaluate at, qprime = q - q_shift is the variable we differentiate against
template <typename Real>
void domegadq_shift(Real &xcomp, Real &ycomp, const SymmetryCoordinates::CartesianPoint2<Real> &q, const SymmetryCoordinates::CartesianPoint2<Real> &q_shift, fadbad::F<Real> (*fadfunc)(fadbad::F<Real>, fadbad::F<Real>))
{
    fadbad::F<double> sr,st;
    auto qdiffsym = SymmetryCoordinates::fromCartesian(q - q_shift);
    sr = qdiffsym.r;
    st = qdiffsym.t;
    sr.diff(0,2);
    st.diff(1,2);
    auto qsym = SymmetryCoordinates::fromCartesian(q);
    sr = qsym.r;
    st = qsym.t;

    auto w = LATTICE_A*1e12*fadfunc(sr, st);
    const double dwdrg = w.d(0);
    const double dwdthetag = w.d(1);
    const double polar_r = sqrt(q.x*q.x + q.y*q.y);
    const double polar_t = atan2(q.y, q.x);
    const double thetag = qsym.t;
    const double modk = polar_r;
    const double dwdr = dwdrg * cos(thetag);
    const double dwdt = dwdrg * (-modk * sin(thetag)) + dwdthetag;
    const double dwdx = dwdr * cos(polar_t) - polar_r * dwdt * sin(polar_t);
    const double dwdy = dwdr * sin(polar_t) + polar_r * dwdt * cos(polar_t);
    xcomp = dwdx;
    ycomp = dwdy;
    return;
}

template <bool splitting, typename Real>
std::vector<SymmetryCoordinates::CartesianPoint2<Real>> CalculateQDoublePrime(const SymmetryCoordinates::CartesianPoint2<Real> &q, const std::vector<SymmetryCoordinates::CartesianPoint2<Real>> &qprime)
{
    std::vector<SymmetryCoordinates::CartesianPoint2<Real>> inp(qprime);

    for(auto &cp : inp) {
        auto const diff = q - cp;
        auto const sp = SymmetryCoordinates::fromCartesian(diff);

        /* inp contains q'' in cartesian coordinates. */
        cp.x = sp.r * std::tan(sp.t) / SymmetryCoordinates::L;
        cp.y = std::acos(std::cos(M_PI*sp.r)) / (M_PI * SymmetryCoordinates::L);
    }

    return inp;
}

template <bool splitting, typename Real>
std::vector<Real> IntegrateTauInverse(const SymmetryCoordinates::CartesianPoint2<Real> &q, const std::vector<SymmetryCoordinates::CartesianPoint2<Real>> &coord, std::vector<Real> &velocity_branch)
{
    auto curve = CalculateQPrime<splitting, Real>(q, coord);
    auto inp = CalculateQDoublePrime<splitting, Real>(q, curve);

    const double hexagon_area = (3.0 * sqrt(3.0) / 2.0) * (2.0 / 3.0) * (2.0 / 3.0);
    const double dqA = hexagon_area / coord.size();

    auto qsym = SymmetryCoordinates::fromCartesian(q);

    std::vector<SymmetryCoordinates::SymmetryPoint2<Real>> curvesym;
    curvesym.reserve(curve.size());

    std::vector<SymmetryCoordinates::SymmetryPoint2<Real>> inpsym;
    inpsym.reserve(curve.size());

    for (size_t i = 0; i < curve.size(); i++) {
        /*
            qsym, curvesym, and inpsym are symmetry coordinates of
            q, q', and q'' respectively.
        */
        curvesym.push_back(fromCartesian(curve[i]));
        inpsym.push_back(fromCartesian(inp[i]));
    }

    Real (*w_functions[])(Real, Real) = {
        w1, w2, w3, w4, w5, w6
    };
    const size_t num_w = sizeof(w_functions)/sizeof(w_functions[0]);
    fadbad::F<double> (*w_fad_functions[])(fadbad::F<double>, fadbad::F<double>) = {
         w1, w2, w3, w4, w5, w6
     };

    // Initialize sums with zero
    std::vector<double> tau_inv_branch(num_w, 0);

    // Large tolerance because omega is in the terahertz range.
    const double reltol = 1e-1;
    const double tol = 1e12*reltol;

    // For each point which has a q' and q'', calculate omega, omega',
    // and omega''. If the combination for a splitting or combining process is
    // within the tolerance, we consider that an intersection and add the
    // contribution to tau inverse.
    for (size_t i = 0; i < curve.size(); i++) {
        for (size_t j = 0; j < num_w; j++) {
            const Real w = 1e12 * w_functions[j](qsym.r, qsym.t);
            for (size_t k = 0; k < num_w; k++) {
                const Real w_prime = 1e12 * w_functions[k](curvesym[i].r, curvesym[i].t);
                double dwprimedqprimex = 0;
                double dwprimedqprimey = 0;
                domegadq(dwprimedqprimex, dwprimedqprimey, curve[i], curvesym[i], w_fad_functions[k]);
                for (size_t l = 0; l < num_w; l++) {
                    const Real w_doubleprime = 1e12 * w_functions[l](inpsym[i].r, inpsym[i].t);
                    double dwppdqprimex = 0;
                    double dwppdqprimey = 0;
                    domegadq_shift(dwppdqprimex, dwppdqprimey, inp[i], inp[i]-curve[i], w_fad_functions[l]);
                    const double vgx = dwprimedqprimex - dwppdqprimex;
                    const double vgy = dwprimedqprimey - dwppdqprimey;
                    const double group_velocity = std::sqrt(vgx*vgx + vgy*vgy);
                    Real delta_w = 0;
                    if (splitting) {
                        delta_w = w - w_prime - w_doubleprime;
                    } else {
                        delta_w = w + w_prime - w_doubleprime;
                    }
                    if (fabs(delta_w) < tol) {
                        Real N = 0;
                        if (splitting) {
                            N = NSplitting(w_prime, w_doubleprime);
                        } else {
                            N = NCombining(w_prime, w_doubleprime);
                        }
                        tau_inv_branch[j] += w*w_prime*w_doubleprime*dqA*N / group_velocity;

                        if (N < 0) {
                            std::cout << "splitting: " << splitting << '\n';
                            std::cout << "w,w',w'': " << w << ',' << w_prime << ',' << w_doubleprime << '\n';
                            std::cout << "N < 0: " << N << ", jkl = " << j << ',' << k << ',' << l << '\n'; 
                        }
                    }
                }
            }
        }
    }

    // We have to scale everything by these factors to get the units correct.
    fadbad::F<double> qx, qy;
    qx = q.x;
    qy = q.y;
    qx.diff(0,2);
    qy.diff(1,2);
    SymmetryCoordinates::CartesianPoint2<fadbad::F<double>> cq= {qx, qy};
    auto sq = SymmetryCoordinates::fromCartesian(cq);
   
    for (size_t i = 0; i < tau_inv_branch.size(); i++) { 
        /*
        fadbad::F<double> sr,st;
        sr = qsym.r;
        st = qsym.t;
        sr.diff(0,2);
        st.diff(1,2);

        auto w = LATTICE_A*1e12*w_fad_functions[i](sr, st);
        const double dwdrg = w.d(0);
        const double dwdthetag = w.d(1);
        double polar_r = sqrt(q.x*q.x + q.y*q.y);
        double polar_t = atan2(q.y, q.x);
        const double rg = qsym.r;
        const double thetag = qsym.t;
        double modk = polar_r;
        double dwdr = dwdrg * cos(thetag);
        double dwdt = dwdrg * (-modk * sin(thetag)) + dwdthetag;
        double dwdx = dwdr * cos(polar_t) - polar_r * dwdt * sin(polar_t);
        double dwdy = dwdr * sin(polar_t) + polar_r * dwdt * cos(polar_t);
        */

        double dwdx = 0;
        double dwdy = 0;
        domegadq(dwdx, dwdy, q, qsym, w_fad_functions[i]);

        auto &tau_inv = tau_inv_branch[i];
        const double velocity = SymmetryCoordinates::L * std::sqrt(dwdx*dwdx + dwdy*dwdy);
        // const double velocity = SymmetryCoordinates::L * std::sqrt(dwdrg*dwdrg + (dwdthetag*dwdthetag) / (modk*modk) - 2*sin(thetag)*dwdrg*dwdthetag);
        // std::cout << "dwdq: (" << dwdx << ", " << dwdy << ")\n";
        // std::cout << "velocity: " << velocity << '\n';
        velocity_branch[i] = velocity;
        const double f = (2 * GRUNEISEN_PARAMETER * GRUNEISEN_PARAMETER * DIRAC_CONSTANT) / (3 * M_PI * GRAPHENE_DENSITY * velocity * velocity);
        const double qscale = 2 * M_PI / LATTICE_A;
        tau_inv *= f * qscale * qscale * GRAPHENE_SAMPLE_LENGTH * GRAPHENE_SAMPLE_LENGTH;
    }

    for (size_t i = 0; i < tau_inv_branch.size(); i++) {
        const double velocity = velocity_branch[i];
        const double specularity = 0.99;
        const double tau_boundary = (velocity / GRAPHENE_SAMPLE_LENGTH) * (1.0 - specularity) / (1.0 + specularity);
        tau_inv_branch[i] += tau_boundary;
    }

    return tau_inv_branch;
}


int main(int argc, char **argv)
{
    if (argc < 5) {
        std::cout << argv[0] << ": GRIDFILE QX QY OUTPUT-PREFIX.\n";
        return 0;
    }
    std::string infile(argv[1]);
    SymmetryCoordinates::CartesianPoint2<double> q = {std::stod(argv[2]), std::stod(argv[3])};

    // std::cout << "Q point is (" << q << ").\n";

    std::vector<SymmetryCoordinates::CartesianPoint2<double>> coord;
    coord.reserve(POINTS_HINT);
    std::ifstream fin(infile);
    while (true) {
        SymmetryCoordinates::CartesianPoint2<double> cp;
        fin >> cp.x;
        fin >> cp.y;
        if (fin.fail()) {
            break;
        }
        coord.push_back(cp);
    }
    // std::cout << "Read " << coord.size() << " grid points.\n";

    std::vector<double> velocity_branch(6);
    auto tau_inv_branch_splitting = IntegrateTauInverse<true, double>(q, coord, velocity_branch);
    auto tau_inv_branch_combining = IntegrateTauInverse<false, double>(q, coord, velocity_branch);

    if (tau_inv_branch_combining.size() != tau_inv_branch_splitting.size()) {
        throw std::runtime_error("Size mismatch in processes.");
    }

    std::vector<double> tau_inv_branch(tau_inv_branch_splitting.begin(), tau_inv_branch_splitting.end());
    for (size_t i = 0; i < tau_inv_branch_combining.size(); i++) {
        tau_inv_branch[i] += tau_inv_branch_combining[i];
    }

    for (size_t i = 0; i < tau_inv_branch.size(); i++) {
        // std::cout << "tau_inv_branch[" << i << "] = " << tau_inv_branch[i] << '\n';
    }

    std::string filename(argv[4]);
    filename = filename+std::to_string(q.x)+'_'+std::to_string(q.y);
    std::ofstream ofs(filename);
    ofs << q.x << ' ' << q.y;
    for (auto tau_inv: tau_inv_branch) {
        if (std::isinf(tau_inv)) {
            ofs << ' ' << 1e300;
        } else {
            ofs << ' ' << tau_inv;
        }
    }

    for (auto velocity: velocity_branch) {
        ofs << ' ' << velocity;
    }

    return 0;
}
