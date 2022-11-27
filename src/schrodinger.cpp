#include "schrodinger.hpp"
#include <armadillo>
#include <cassert>
#include <random>
#include <iostream>
#include <math.h>

std::complex<double> i_unit = std::complex<double>(0., 1.);

// Constructor
Schrodinger::Schrodinger(const int M_given, const double h_given, const double dt_given)
{

    M = M_given;
    h = h_given;
    dt = dt_given;

    // set r parameter
    r = 0.5 * i_unit * dt / (h * h);
}

// Set the potential matrix
void Schrodinger::set_potential(arma::mat &V_given)
{
    // wall between 0.49 and 0.51
    V = &V_given;
    (*V).fill(0.);


    
    double x, y, V0 = 1e+10;
    // construct vertical walls
    for (int i{}; i < M - 2; i++)
    {
        x = i * h;
        if (x >= 0.49 && x<= 0.51)
        {
            (*V).col(i).fill(V0);   
        }
    }
    // make slits in the wall, just double slit for now...
    for(int j{}; j < M-2; j++)
    {
        y = j * h;
        if(y >= 0.425 && y <= 0.475 )
            (*V).row(j).fill(0.);
        if(y >= 0.525 && y <= 0.575 )
            (*V).row(j).fill(0.);

    }

}

// Set A and B matrices
void Schrodinger::set_AB()
{
    // temp variable
    double M_2 = M - 2;

    arma::cx_vec a = arma::cx_vec((M - 2) * (M - 2)).fill(0.);
    arma::cx_vec b = arma::cx_vec((M - 2) * (M - 2)).fill(0.);
    for (int i{}; i < M_2; i++)
    {
        for (int j{}; j < M_2; j++)
        {
            a(single_index(i, j)) = 1. + 4. * r + i_unit * 0.5 * dt * (*V)(i, j);
            b(single_index(i, j)) = 1. - 4. * r - i_unit * 0.5 * dt * (*V)(i, j);
        }
    }

    // construct matrices
    set_A(r, a);
    set_B(r, b);
}

// Switch from i,j indices to single index
int Schrodinger::single_index(const int i, const int j)
{
    return i + (M - 2) * j;
}

// Set up the A matrix for Crank-Nicholson method
void Schrodinger::set_A(const arma::cx_double r, const arma::cx_vec &a)
{
    int M_2 = M - 2;

    // Start from identity matrix
    A.eye(M_2 * M_2, M_2 * M_2);

    // Loop that fills rows with repating "block"
    for (int i{}; i < M_2 * M_2; i++)
    {

        A(i, i) = a(i);

        if ((i + 1) % M_2 != 0)
        {
            A(i, i + 1) = -r;
            A(i + 1, i) = -r;
        }
    }

    // loop that fills the outer off diagonals -r
    int i_max = M_2 * (M_2 - 1);
    for (int i = 0; i < i_max; i++)
    {
        A(i, i + M_2) = -r;
        A(i + M_2, i) = -r;
    }
}

// Set up the B matrix for Crank-Nicholson method
void Schrodinger::set_B(const arma::cx_double r, const arma::cx_vec &b)
{

    int M_2 = M - 2;

    // Start from identity matrix
    B.eye(M_2 * M_2, M_2 * M_2);

    // Loop that fills rows with repating "block"
    for (int i{}; i < M_2 * M_2; i++)
    {

        B(i, i) = b(i);

        if ((i + 1) % M_2 != 0)
        {
            B(i, i + 1) = r;
            B(i + 1, i) = r;
        }
    }

    // loop that fills the outer off diagonals +r
    int i_max = M_2 * (M_2 - 1);
    for (int i = 0; i < i_max; i++)
    {
        B(i, i + M_2) = r;
        B(i + M_2, i) = r;
    }
}

// Set up the initial state in a gaussian wave packet
void Schrodinger::set_U(const double x_c, const double y_c, const double sgm_x, const double sgm_y, const double p_x, const double p_y, arma::cx_vec &u_given)
{
    u = &u_given;
    double sgm_x2 = sgm_x * sgm_x;
    double sgm_y2 = sgm_y * sgm_y;
    double x = 0, y = 0;

    arma::cx_double normalization = arma::cx_double(0., 0.);
    for (int i = 1; i < M - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            x = i * h;
            y = j * h;
            double delta_x = x - x_c, delta_y = y - y_c;
            (*u)(single_index(i - 1, j - 1)) = std::exp(-1. * delta_x * delta_x / (2 * sgm_x2) - 1. * delta_y * delta_y / (2 * sgm_y2) + i_unit * p_x * delta_x + i_unit * p_y * delta_y);
            normalization += std::abs((*u)(single_index(i - 1, j - 1))); // sum of |u|^2 for normalization
        }
        // Above we can just compute the u vector, and if we want to store to file
        // the full solution matrix U maybe we can just "reshape" the vector back
        // into a matrix with .reshape() from armadillo ? We'll see... tocca.
    }
    // compute renormalization
    normalization /= (M - 2) * (M - 2);
    normalization = std::sqrt(normalization);

    // normalize
    (*u) = (*u) / normalization;
}

// Evolve the system
void Schrodinger::evolve()
{
    arma::cx_vec b = arma::cx_vec(M - 2);

    // calculate b vector
    b = B * (*u);

    // solve matrix equation A*u = b
    (*u) = arma::solve(A, b);
}

// Print on scren and file the system at the current state
void Schrodinger::print_data()
{
    int M_2 = M - 2;
    arma::cx_mat U; // U matrix
    U.set_size(M, M);
    for (int i = 1; i < M - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            U(i, j) = (*u)(single_index(i - 1, j - 1));
        }
    }
    // here in the loop we can directly print to file
    // without the need of having the actual matrix U saved in memory!

    U.col(0).fill(std::complex<double>(0., 0.));
    U.col(M - 1).fill(std::complex<double>(0., 0.));
    U.row(0).fill(std::complex<double>(0., 0.));
    U.row(M - 1).fill(std::complex<double>(0., 0.));

    U.print();
    //U.save("data/state_evolved.txt", arma::raw_ascii);

}