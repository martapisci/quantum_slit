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
    int M = 40;            // number of points along one side
    double h = 1./M;     // spatial step
    double dt = 0.000025; // time step
    double T = 0.008;     // final time
    double x_c = 0.25;    // Coordinate (x) of the center of the wave packet
    double y_c = 0.5;     // Coordinate (y) of the center of the wave packet;
    double sgm_x = 0.05;  // Initial width (x) of the wave packet
    double sgm_y = 0.05;  // Initial width (y) of the wave packet
    double p_x = 200.;    // Initial momentum (x) of the wave packet
    double p_y = 0.;      // Initial momentum (y) of the wave packet
    // ----------------------------------------------------------

    // setup imaginary unit
    // std::complex<double> i_unit = std::complex<double>(0., 1.);

    int M_2 = M - 2;
    // initialize matrices and vectors
    arma::mat V; // Potential matrix
    V.set_size(M_2, M_2);
    arma::cx_vec u = arma::cx_vec(M_2 * M_2);

    // initialize lattice
    Schrodinger lattice = Schrodinger(M, h, dt);
    // set potential
    lattice.set_potential(V);
    // set Crank-Nicholson matrices
    lattice.set_AB();
    // set initial state
    lattice.set_U(x_c, y_c, sgm_x, sgm_y, p_x, p_y, u);

    // print to file the potential
    V.save("data/potential.txt", arma::raw_ascii);

    // evolve!
    lattice.evolve();

    //lattice.print_data();

    // debug line
    // std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << std::endl;

    std::cout << "All good!" << std::endl;
    // all is good
    return 0;
}