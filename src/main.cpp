#include <iostream>
#include <stdlib.h>
#include "utils.hpp"
#include "schrodinger.hpp"
#include <armadillo>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>

int main()
{
    // ----------------------- PARAMETERS -----------------------
    int M = 5;        // number of points along one side
    double h = 0.01;  // spatial step
    double dt = 0.01; // time step
    // ----------------------------------------------------------

    // setup imaginary unit
    // std::complex<double> i_unit = std::complex<double>(0., 1.);

    int M_2 = M - 2;
    // initialize matrices and vectors
    arma::mat V = arma::mat(M_2, M_2).randu(); // Potential matrix
    V = 1000* V;

    // initialize lattice
    Schrodinger lattice = Schrodinger(M, h, dt, V);


    // debug line
    // std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << std::endl;

    std::cout << "All good so far... the class constructor seems to work properly :)" << std::endl;
    // all is good
    return 0;
}