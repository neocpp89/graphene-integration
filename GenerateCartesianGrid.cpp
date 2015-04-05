#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstdlib>

#include "SymmetryCoordinates.hpp"

using Point2 = SymmetryCoordinates::CartesianPoint2<double>;

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

std::vector<Point2> UniformRadialAndTheta(size_t num_radial, size_t num_theta)
{
    std::vector<Point2> output(num_radial*num_theta);

    /*
        Scan radially (outer) and angularly (inner) to create a uniform grid
        in (r,theta). Note that the point density per area is not constant
        when done this way.
    */
    const double rsqrt3 = 1.0 / std::sqrt(3.0);
    for (size_t m = 0; m < num_radial; m++) {
        const double rho = rsqrt3 * m / (num_radial-1);
        for (size_t n = 0; n < num_theta; n++) {
            const double theta = 2*M_PI * n/(num_theta-1);
            output[m*num_theta + n] = F(rho, theta);
        }
    }

    return output;
}

std::vector<Point2> UniformCartesian(size_t Nx, size_t Ny)
{
    std::vector<Point2> output;
    output.reserve(Nx*Ny);

    const double hx = 2.0 / (Nx-1);
    const double hy = 2.0 / (Ny-1);

    // Generate a cartesian grid a single point at a time
    for (size_t i = 0; i < Nx; i++) {
        for (size_t j = 0; j < Ny; j++) {
            Point2 p = {i*hx - 1.0, j*hy - 1.0};
            auto const sp = SymmetryCoordinates::fromCartesian(p);

            /*
                Check if the point is within the Brillouin zone by using the
                symmetry coordinates.
            */
            if (sp.r <= 1) {
                output.push_back(p);
            }
        }
    }

    return output;
}

int main(int argc, char **argv)
{
    std::string outfile;
    bool CartesianMode = true;
    if (argc < 4) {
        std::cout << argv[0] << ": OUTPUTFILE N M [type].\n";
        std::cout << "    OUTPUTFILE - where to write points within Brillouin zone.\n";
        std::cout << "    N - number of x points or radial points.\n";
        std::cout << "    M - number of y points or theta points.\n";
        std::cout << "    type - if 'R', use radial mode, otherwise use cartesian mode (default).\n";
        return 0;
    } else {
        if (argc >= 5) {
            if (std::string(argv[4]) == std::string("R")) {
                CartesianMode = false;
            } else {
                std::cout << "Unknown mode (leave blank for Cartesian).\n";
            }
        }
    }
    outfile = argv[1];
    const size_t N = std::stoul(argv[2]);
    const size_t M = std::stoul(argv[3]);

    std::vector<Point2> output;

    if (CartesianMode) {
        output = UniformCartesian(N, M);
    } else {
        output = UniformRadialAndTheta(N, M);
    }

    std::ofstream fout(outfile);
    for (auto const &item : output) {
        fout << item.x << ' '<< item.y << '\n';
    }
    
    return 0;
}
