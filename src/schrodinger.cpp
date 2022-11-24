#include "schrodinger.hpp"
#include <armadillo>
#include <cassert>
#include <random>
#include <iostream>

// Constructor
Schrodinger::Schrodinger(const int M_given, const double h_given, const double dt_given, arma::mat &V_given)
{

    M = M_given;
    h = h_given;
    dt = dt_given;

    // temp variable
    double M_2 = M - 2;

    // setup imaginary unit
    std::complex<double> i_unit = std::complex<double>(0., 1.);
    
    // set r parameter
    r = 0.5 * i_unit * dt / (h * h);

    arma::cx_vec a = arma::cx_vec((M - 2) * (M - 2)).fill(0.);
    arma::cx_vec b = arma::cx_vec((M - 2) * (M - 2)).fill(0.);
    for (int i{}; i < M_2; i++)
    {
        for (int j{}; j < M_2; j++)
        {
            a(single_index(i, j)) = 1. + 4. * r + i_unit * 0.5 * dt * V_given(i, j);
            b(single_index(i, j)) = 1. - 4. * r - i_unit * 0.5 * dt * V_given(i, j);
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

    // loop that fills the outer off diagonals -r
    int i_max = M_2 * (M_2 - 1);
    for (int i = 0; i < i_max; i++)
    {
        B(i, i + M_2) = r;
        B(i + M_2, i) = r;
    }
}

// Set up the initial state in a gaussian wave packet
void Schrodinger::set_U(const double x_c, const double y_c, const double sgm_x, const double sgm_y, const double p_x, const double p_y)
{
}

// Evolve the system by a spatial step h and a time step dt
void Schrodinger::evolve(const double h, const double dt)
{
}
