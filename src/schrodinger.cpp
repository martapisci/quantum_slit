#include "schrodinger.hpp"
#include "utils.hpp"
#include <armadillo>
#include <cassert>
#include <random>
#include <iostream>
#include <math.h>
#include <complex>

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
void Schrodinger::set_potential(arma::sp_mat &V_given, const int &switch_given)
{
    // wall between 0.49 and 0.51
    V = &V_given;

    double x, y, V0 = 1e+10;

    if (switch_given == 0)
        V0 = 0.;

    // construct vertical walls
    for (int i{}; i < M - 2; i++)
    {
        x = i * h;
        if (x >= 0.49 && x <= 0.51)
        {
            (*V).row(i).fill(V0);
        }
    }
    // make slits in the wall, just double slit for now...
    for (int j{}; j < M - 2; j++)
    {
        y = j * h;
        if (y >= 0.425 && y <= 0.475)
            (*V).col(j).fill(0.);
        if (y >= 0.525 && y <= 0.575)
            (*V).col(j).fill(0.);
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
            std::complex<double> k = std::complex<double>((*V)(i, j), 0.);
            a(single_index(i, j)) = 1. + 4. * r + i_unit * 0.5 * dt * k;
            b(single_index(i, j)) = 1. - 4. * r - i_unit * 0.5 * dt * k;
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
    // For every submatrix...
    for(int i_block{}; i_block< M_2; i_block++)
    {
        int i_tmp = M_2*i_block;
        for(int i{}; i < M_2; i++)
        {
            int k = i_tmp + i;
            // fill the submatrix with diagonals
            A(k, k) = a(k); 
            // and fill with off diagonal r's except the last step
            if(i != M_2 -1)
            {
                A(k, k + 1) = -r;
                A(k + 1, k) = -r;
            }
        }
    }
    // Now fill outer off-diagonals
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
    // For every submatrix...
    for(int i_block{}; i_block< M_2; i_block++)
    {
        int i_tmp = M_2*i_block;
        for(int i{}; i < M_2; i++)
        {
            int k = i_tmp + i;
            // fill the submatrix with diagonals
            B(k, k) = b(k); 
            // and fill with off diagonal r's except the last step
            if(i != M_2 -1)
            {
                B(k, k + 1) = r;
                B(k + 1, k) = r;
            }
        }
    }
    // Now fill outer off-diagonals
    int i_max = M_2 * (M_2 - 1);
    for (int i = 0; i < i_max; i++)
    {
        B(i, i + M_2) = r;
        B(i + M_2, i) = r;
    }
}

// Set up the initial state in a gaussian wave packet
void Schrodinger::set_initial_state(const double x_c, const double y_c, const double sgm_x, const double sgm_y, const double p_x, const double p_y, arma::cx_vec &u_given)
{
    u = &u_given;
    double sgm_x2 = sgm_x * sgm_x;
    double sgm_y2 = sgm_y * sgm_y;
    double x = 0, y = 0;

    double normalization = 0.;
    for (int i = 1; i < M - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            x = i * h;
            y = j * h;
            double delta_x = x - x_c, delta_y = y - y_c;
            std::complex<double> tmp = std::exp(-1. * delta_x * delta_x / (2. * sgm_x2) - 1. * delta_y * delta_y / (2. * sgm_y2) + i_unit * p_x * delta_x + i_unit * p_y * delta_y);
            (*u)(single_index(i - 1, j - 1)) = tmp;
            normalization += std::real(tmp*std::conj(tmp));

        }
    }

    // normalize
    (*u) = (*u) /(std::sqrt(normalization));
}

// Evolve the system
void Schrodinger::evolve()
{
    arma::cx_vec b = arma::cx_vec(M - 2);

    // calculate b vector
    b = B * (*u);

    // solve matrix equation A*u = b
    arma::spsolve( (*u), A, b);
}

//Update the U state matrix
void Schrodinger::U(arma::cx_mat &U)
{
    // copy state
    arma::cx_mat u_tmp = *u;

    // shape it as a matrix
    u_tmp.reshape(M-2,M-2);
    
    // update the 'inner' state matrix
    U(arma::span(1, M - 2), arma::span(1, M - 2)) = u_tmp;
}
// Calculate current total probability
void Schrodinger::probability(double &prob)
{
    double probability = 0.;
    for (int i = 1; i < M - 1; i++)
    {
        for (int j = 1; j < M - 1; j++)
        {
            std::complex<double> tmp = (*u)(single_index(i - 1, j - 1));
            probability += std::real(tmp*std::conj(tmp)); // sum of |u_ij|^2
        }
    }

    prob = probability;
}