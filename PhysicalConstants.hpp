#ifndef PHYSICAL_CONSTANTS_HPP
#define PHYSICAL_CONSTANTS_HPP
#include <cmath>

const double PLANCK_CONSTANT = 6.626e-34; // m^2 kg / s
const double DIRAC_CONSTANT = PLANCK_CONSTANT / (2 * M_PI); // m^2 kg / s
const double BOLTZMANN_CONSTANT = 1.3806e-23; // m^2 kg s^-2 K^-1
const double GRAPHENE_DENSITY = 7.6e-7; // areal density kg / m^2
const double CARBON_MASS = 1.994e-23; // kg
const double LATTICE_A = 2.46e-10; // m
const double GRUNEISEN_PARAMETER = 1.5; 
const double GRAPHENE_THICKNESS = 3.4e-10; // m
const double GRAPHENE_SAMPLE_LENGTH = 5e-6; // m

const size_t NUMBER_OF_BRANCHES = 6;

#endif //PHYSICAL_CONSTANTS_HPP
