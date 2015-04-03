#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>

#include "fadiff.h"
#include "SymmetryCoordinates.hpp"
#include "Interpolations.hpp"

const size_t POINTS_HINT = 4000000;

double delta(double arg, double sigma=10.0)
{
    return std::exp(arg*arg/(2*sigma*sigma));
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cout << "Needs cartesian grid file and output file.\n";
        return 0;
    }

    std::vector<SymmetryCoordinates::CartesianPoint2> coord;
    coord.reserve(POINTS_HINT);
    std::ifstream fin(argv[1]);
    while (true) {
        SymmetryCoordinates::CartesianPoint2 cp;
        fin >> cp.x;
        fin >> cp.y;
        if (fin.fail()) {
            break;
        }
        coord.push_back(cp);
    }

    /*
    for (auto const &cp : coord) {
        std::cout << cp.x << ' ' << cp.y << '\n';
    }
    */

    SymmetryCoordinates::CartesianPoint2 q = {0.15,0.15};
   
    std::vector<SymmetryCoordinates::CartesianPoint2> curve;
    curve.reserve(coord.size());

    for (auto const &cp : coord) {
        SymmetryCoordinates::CartesianPoint2 diff = q - cp;
        SymmetryCoordinates::SymmetryPoint2 sp = SymmetryCoordinates::fromCartesian(diff);
        if (sp.r > 1) { 
            /* curve contains q' in cartesian coordinates. */
            curve.push_back(cp);
        }
    }

    std::vector<SymmetryCoordinates::CartesianPoint2> inp(curve);

    for(auto &cp : inp) {
        auto diff = q - cp;
        auto sp = SymmetryCoordinates::fromCartesian(diff);

        /* inp contains q'' in cartesian coordinates. */
        cp.x = sp.r * std::tan(sp.t) / SymmetryCoordinates::L;
        cp.y = std::acos(std::cos(M_PI*sp.r)) / (M_PI * SymmetryCoordinates::L);
    }

    auto qsym = SymmetryCoordinates::fromCartesian(q);

    std::vector<SymmetryCoordinates::SymmetryPoint2> curvesym;
    curvesym.reserve(curve.size());

    std::vector<SymmetryCoordinates::SymmetryPoint2> inpsym;
    inpsym.reserve(curve.size());

    for (size_t i = 0; i < curve.size(); i++) {
        /*
            qsym, curvesym, and inpsym are symmetry coordinates of
            q, q', and q'' respectively.
        */
        curvesym.push_back(fromCartesian(curve[i]));
        inpsym.push_back(fromCartesian(inp[i]));
    }

    double (*w_functions[])(double, double) = {
        w1, w2, w3, w4, w5, w6
    };
    fadbad::F<double> (*w_fad_functions[])(fadbad::F<double>, fadbad::F<double>) = {
        w1, w2, w3, w4, w5, w6
    };
    const size_t num_w = sizeof(w_functions)/sizeof(w_functions[0]);

    fadbad::F<double> qr, qt;
    qr = qsym.r;
    qr.diff(0,2);
    qt = qsym.t;
    qt.diff(1,2);

    for (size_t j = 0; j < num_w; j++) {
    for (size_t k = 0; k < num_w; k++) {
    for (size_t l = 0; l < num_w; l++) {
        std::vector<SymmetryCoordinates::CartesianPoint2> scurve;
        std::vector<SymmetryCoordinates::CartesianPoint2> gradient;
        scurve.reserve(curve.size());
        const double tol = 1e-1;
        for (size_t i = 0; i < curve.size(); i++) {
            double val = w_functions[j](qsym.r, qsym.t)
                         - w_functions[k](curvesym[i].r, curvesym[i].t)
                         - w_functions[l](inpsym[i].r, inpsym[i].t);
            if (std::fabs(val) < tol) {
                fadbad::F<double> cr, ct;
                fadbad::F<double> ir, it;
                cr = curvesym[i].r;
                cr.diff(0,2);
                ct = curvesym[i].t;
                ct.diff(1,2);
                ir = inpsym[i].r;
                ir.diff(0,2);
                it = inpsym[i].t;
                it.diff(1,2);

                auto f = w_fad_functions[j](qr, qt);
                auto g = -w_fad_functions[k](cr, ct);
                auto h = -w_fad_functions[l](ir, it);
                double dfdx = f.d(0) * std::cos(qsym.t) - f.d(0) * qsym.r * std::sin(qsym.t);
                double dfdy = f.d(1);
                double dhdx = h.d(0) * std::cos(inpsym[i].t) - h.d(0) * inpsym[i].r * std::sin(inpsym[i].t);
                double dhdy = h.d(1);
                double dgdx = g.d(0) * std::cos(curvesym[i].t) - g.d(0) * curvesym[i].r * std::sin(curvesym[i].t);
                double dgdy = g.d(1);
                scurve.push_back(curve[i]);
                gradient.push_back({dfdx+dgdx+dhdx, dfdy+dgdx+dhdx});
                std::cout << dfdx << ", " << dfdy << '\n';
            }
        }

        if (scurve.size() > 0) {
            std::string filename("output/scurve");
            std::ofstream foutsc(filename+std::to_string(j)+std::to_string(k)+std::to_string(l));
            for (auto const &cp : scurve) {
                foutsc << cp.x << ' ' << cp.y << '\n';
            }
            std::ofstream foutgrad(filename+"grad"+std::to_string(j)+std::to_string(k)+std::to_string(l));
            for (auto const &grad : gradient) {
                const double t = grad.y;
                const double r = grad.x;
                const double x = r*std::cos(t);
                const double y = r*std::sin(t);
                foutgrad << x << ' ' << y << '\n';
            }
        }
        std::cout << "Finished: (" << j << ", " << k << ", " << l << ")\n";
    }
    }
    }

    std::ofstream fout(argv[2]);
    for (auto const &cp : curve) {
        fout << cp.x << ' ' << cp.y << '\n';
    }

    std::ofstream fouti("inp");
    for (auto const &cp : inp) {
        fouti << cp.x << ' ' << cp.y << '\n';
    }


    return 0;
}
