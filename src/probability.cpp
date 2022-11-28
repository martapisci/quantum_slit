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
    int M;               // number of points along one side
    std::cout << "Insert M: " << std::endl;
    std::cin >> M;
    double h = 1./M;        // spatial step
    double dt = 0.000025;     // time step
    double T = 0.008;         // final time
    double x_c = 0.25;        // Coordinate (x) of the center of the wave packet
    double y_c = 0.5;         // Coordinate (y) of the center of the wave packet;
    double sgm_x = 0.05;      // Initial width (x) of the wave packet
    double sgm_y = 0.05;      // Initial width (y) of the wave packet
    double p_x = 200.;        // Initial momentum (x) of the wave packet
    double p_y = 0.;          // Initial momentum (y) of the wave packet
    int switch_potential; // Potential switch
    std::cout << "Insert 0 for no potential or else for double slit: " << std::endl;
    std::cin >> switch_potential;

    // ----------------------------------------------------------

    // total current probabilty
    double probability_now;
    // Print parameters
    int width = 10;
    int prec = 10;

    // setup imaginary unit
    // std::complex<double> i_unit = std::complex<double>(0., 1.);

    // handy temp. variable
    int M_2 = M - 2;

    // initialize matrices and vectors
    arma::sp_mat V = arma::sp_mat(M_2, M_2); // Potential matrix
    // V.set_size(M_2, M_2);
    arma::cx_vec u = arma::cx_vec(M_2 * M_2);
    // initialize lattice
    Schrodinger lattice = Schrodinger(M, h, dt);

    // set potential
    lattice.set_potential(V, switch_potential);

    // set Crank-Nicholson matrices
    lattice.set_AB();

    // set initial state
    lattice.set_initial_state(x_c, y_c, sgm_x, sgm_y, p_x, p_y, u);

    // check that probability is 1
    lattice.probability(probability_now);
    std::cout << "Total initial probability is " << scientific_format(probability_now, width, prec) << std::endl;

    // Set file for probability
    std::string potential_str = switch_potential ? "_slit" : "_no_slit";
    std::string filename = "data/probability_M" + std::to_string(M) + potential_str+ ".txt";
    std::ofstream ofile;
    ofile.open(filename);
    // Set file for initial state
    std::string filename2 = "data/state0_M" + std::to_string(M) + ".txt";
    std::ofstream ofile2;
    ofile2.open(filename2);

    // Print to file the initial state
    lattice.print_data(ofile2, width, prec);

    // time loop
    int T_steps = T / dt;
    double t = 0.;
    ofile << scientific_format(t, width, prec) << " " << scientific_format(probability_now - 1., width, prec) << std::endl;
    for (int i = 1; i < T_steps; i++)
    {
        t = i * dt;
        std::cout << "Now calculating t=" << t << std::endl;
        lattice.evolve();
        lattice.probability(probability_now);
        ofile << scientific_format(t, width, prec) << " " << scientific_format(probability_now - 1., width, prec) << std::endl;
    }

    // Print to file the evolved state
    // lattice.print_data(ofile, width, prec);

    // Close file
    ofile.close();
    ofile2.close();

    // debug line
    // std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << std::endl;

    std::cout << "All good!" << std::endl;
    // all is good
    return 0;
}