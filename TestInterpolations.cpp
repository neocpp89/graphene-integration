#include <cmath>
#include <iostream>
#include <fstream>
#include "Interpolations.hpp"

int main()
{
    const double r_max = 1;
    const double t_max = M_PI/6.0;
    const size_t N = 40;
    std::ofstream f1("w1out");
    std::ofstream f2("w2out");
    std::ofstream f3("w3out");
    std::ofstream f4("w4out");
    std::ofstream f5("w5out");
    std::ofstream f6("w6out");
    for (size_t i = 0; i < N; i++) {
        const double r = r_max*i/(N-1);
        for (size_t j = 0; j < N; j++) {
            const double t = t_max*j/(N-1);
            f1 << r*std::cos(t) << ' ' << r*std::sin(t) << ' ' << w1(r, t) << '\n';
            f2 << r*std::cos(t) << ' ' << r*std::sin(t) << ' ' << w2(r, t) << '\n';
            f3 << r*std::cos(t) << ' ' << r*std::sin(t) << ' ' << w3(r, t) << '\n';
            f4 << r*std::cos(t) << ' ' << r*std::sin(t) << ' ' << w4(r, t) << '\n';
            f5 << r*std::cos(t) << ' ' << r*std::sin(t) << ' ' << w5(r, t) << '\n';
            f6 << r*std::cos(t) << ' ' << r*std::sin(t) << ' ' << w6(r, t) << '\n';
        }
    }
    

    return 0;
}
