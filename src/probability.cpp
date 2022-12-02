#include <iostream>
#include <stdlib.h>
#include "utils.hpp"
#include "schrodinger.hpp"
#include <armadillo>
#include <algorithm>
#include <string>
#include <fstream>
#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>

int main()
{
    // Set the filename
    std::string filename_input = "../src/parameters.txt";

    // A vector of vectors to store the rows in the input file
    std::vector<double> input_data;

    // Create a filestream instance "myfile" and use it to read the file
    std::fstream parameters;
    parameters.open(filename_input);
    if (parameters.is_open()) // This checks that the file was opened OK
    {
        // Some temporary variables we'll use
        std::string line;
        double tmp;
        double M_, dt_, T_, x_c_, y_c_, sgm_x_, sgm_y_, p_x_, p_y_;
        int switch_potential_;
        // Read file line by line
        while (std::getline(parameters, line))
        {
            // Skip lines with "#" at the first position
            if (line.at(0) == '#')
            {
                continue;
            }
            else
            {

                std::stringstream mysstream(line);
                mysstream >> tmp;
                input_data.push_back(tmp);
            }
        }
    }
    else
    {
        std::cout << "Unable to open the file " << filename_input;
    }

    // Close the input file
    parameters.close();
    // ----------------------- PARAMETERS -----------------------
    int M = input_data.at(0);                     // number of points along one side
    double h = 1. / M;                            // spatial step
    double dt = input_data.at(1);                 // time step
    double T = input_data.at(2);                  // final time
    double x_c = input_data.at(3);                // Coordinate (x) of the center of the wave packet
    double y_c = input_data.at(4);                // Coordinate (y) of the center of the wave packet;
    double sgm_x = input_data.at(5);              // Initial width (x) of the wave packet
    double sgm_y = input_data.at(6);              // Initial width (y) of the wave packet
    double p_x = input_data.at(7);                // Initial momentum (x) of the wave packet
    double p_y = input_data.at(8);                // Initial momentum (y) of the wave packet
    int switch_potential = (int)input_data.at(9); // Potential switch

    // ----------------------------------------------------------

    // total current probabilty
    double probability_now;
    // Print parameters
    int width = 10;
    int prec = 10;

    // handy temp. variable
    int M_2 = M - 2;

    // initialize matrices and vectors
    arma::cx_mat U = arma::cx_mat(M, M, arma::fill::zeros); // state matrix
    arma::sp_mat V = arma::sp_mat(M_2, M_2);                // Potential matrix
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
    std::string filename = "data/probability_M" + std::to_string(M) + potential_str + ".txt";
    std::ofstream ofile;
    ofile.open(filename);

    // time loop
    int T_steps = T / dt;
    // store data in cubes
    arma::cube U_2 = arma::cube(M, M, T_steps+1, arma::fill::zeros);
    arma::cube U_re = arma::cube(M, M, T_steps+1, arma::fill::zeros);
    arma::cube U_im = arma::cube(M, M, T_steps+1, arma::fill::zeros);
    // initial U
    lattice.U(U);

    U_2.slice(0) = arma::real(U % arma::conj(U));
    U_re.slice(0) = arma::real(U);
    U_im.slice(0) = arma::imag(U);
    double t = 0.;
    ofile << scientific_format(t, width, prec) << " " << scientific_format(probability_now - 1., width, prec) << std::endl;
    for (int i = 1; i <= T_steps; i++)
    {
        t = i * dt;
        std::cout << "Now calculating t=" << t << std::endl;
        lattice.evolve();
        lattice.probability(probability_now);
        ofile << scientific_format(t, width, prec) << " " << scientific_format(probability_now - 1., width, prec) << std::endl;
        // update U
        lattice.U(U);
        U_2.slice(i) = arma::real(U % arma::conj(U));
        U_re.slice(i) = arma::real(U);
        U_im.slice(i) = arma::imag(U);
    }

    // Close file
    ofile.close();

    // save state files
    U_2.save("./data/modulus.bin", arma::arma_binary);
    U_re.save("./data/re.bin", arma::arma_binary);
    U_im.save("./data/im.bin", arma::arma_binary);

    // debug line
    // std::cout << "DEBUG: " << __FILE__ << " " << __LINE__ << std::endl;

    std::cout << "All good!" << std::endl;
    // all is good
    return 0;
}