#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>

#include "SymmetryCoordinates.hpp"

using Point2 = SymmetryCoordinates::CartesianPoint2;

Point2 reflect(const Point2 &p, double theta)
{
    const double ct = std::cos(theta);
    const double st = std::sin(theta);
    const double x = 2*p.y*ct*st + p.x*(ct*ct - st*st);
    const double y = 2*p.x*ct*st + p.y*(-ct*ct + st*st);
    Point2 r = {x, y};
    return r;
}

Point2 F(double rho, double theta, const double n = 6)
{
    Point2 p;
    double t = std::acos(std::cos(n*(theta + M_PI_2)))/n;
    p.x = rho * std::cos(theta) / std::cos(t);
    p.y = rho * std::sin(theta) / std::cos(t);
    return p;
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        std::cout << "Needs number of radial points and file to output.\n";
        return 0;
    }

    

    const size_t N = std::atoi(argv[1]);
    const size_t M = N;

    std::vector<Point2> output;
    output.reserve(M*N);

    const double rsqrt3 = 1.0 / std::sqrt(3.0);

    for (size_t m = 0; m < M; m++) {
        const double rho = rsqrt3 * m / (M-1);
        // const size_t num_theta_points = N*m/(M-1)+10;
        const size_t num_theta_points = N;
        for (size_t n = 0; n < num_theta_points; n++) {
            const double theta = 2*M_PI * n/(num_theta_points-1);
            // const double theta_off = 2*M_PI*static_cast<double>(std::rand()) / RAND_MAX;
            output.push_back(F(rho, theta));
        }
    }
#if 0

    for (size_t n = 0; n < N; n++) {
        for (size_t m = 0; m < M; m++) {
            Point2 p;
            p.x = (m)*(n)/(3.0*(M-1)*(N-1));
            p.y = rsqrt3*m/(M-1);
            output[n*M + m] = p;
        }
    }

    /*
    std::ofstream fout(argv[argc-1]);
    for (size_t n = 0; n < N; n++) {
        for (size_t m = 0; m < M; m++) {
            fout << output[n*M + m].x << ' '<< output[n*M + m].y << '\n';
        }
    }
    */

    std::vector<Point2> fulloutput = output;
    fulloutput.reserve(12*M*N);

    // mirror original data across the y axis
    for (size_t n = 1; n < N; n++) {
        for (size_t m = 0; m < M; m++) {
            Point2 p;
            p.x = -output[n*M + m].x;
            p.y = output[n*M + m].y;
            fulloutput.push_back(p);
        }
    }

    // rotate wedge through all 6 angles
    for (size_t n = 1; n < N; n++) {
        for (size_t m = 0; m < M; m++) {
            Point2 p;
            p.x = -output[n*M + m].x;
            p.y = output[n*M + m].y;
            double t = std::atan2(p.y, p.x);
            fulloutput.push_back(p);
        }
    }

    std::ofstream fout(argv[argc-1]);
    /*
    for (size_t i = 0; i < fulloutput.size(); i++) {
        fout << fulloutput[i].x << ' '<< fulloutput[i].y << '\n';
    }
    */

    std::set<Point2> s(fulloutput.begin(), fulloutput.end());

    for (auto const &item : s) {
        fout << item.x << ' '<< item.y << '\n';
    }
#endif
    std::ofstream fout(argv[argc-1]);
    for (auto const &item : output) {
        fout << item.x << ' '<< item.y << '\n';
    }
    
    return 0;
}
