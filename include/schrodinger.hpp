// include guard
#ifndef __schrodinger_hpp__
#define __schrodinger_hpp__

#include <armadillo>
#include <stdio.h>

class Schrodinger
{
private:
    double h;               // spatial step
    double dt;              // time step
    std::complex<double> r; // r parameter
    arma::cx_mat A;         // A matrix
    arma::cx_mat B;         // B matrix
    // Maybe convenient option to make matrices sp_cx_mat from Armadillo. Look documentation...

public:
    int M;           // number of points along one side
    arma::cx_vec *u; // vectorized matrix with current state
    // arma::Mat<std::complex<double>> *U; // U matrix

    // Constructors
    Schrodinger(const int M_given, const double h, const double dt, arma::mat &V_given);

    /**
     * @brief Switch from i,j indices to single index
     * @param i Row index
     * @param j Column index
     */
    int single_index(const int i, const int j);

    /**
     * @brief Set up the A matrix for Crank-Nicholson method
     * @param r Complex number, off-diagonal element
     * @param a Vector full of diagonal elements
     */
    void set_A(const arma::cx_double r, const arma::cx_vec &a);

    /**
     * @brief Set up the B matrix for Crank-Nicholson method
     * @param r Complex number, off-diagonal element
     * @param b Vector full of diagonal elements
     */
    void set_B(const arma::cx_double r, const arma::cx_vec &b);

    /**
     * @brief Set up the initial state in a gaussian wave packet
     * @param x_c Coordinate (x) of the center of the wave packet
     * @param y_c Coordinate (y) of the center of the wave packet
     * @param sgm_x Initial width (x) of the wave packet
     * @param sgm_y Initial width (y) of the wave packet
     * @param p_x Initial momentum (x) of the wave packet
     * @param p_y Initial momentum (y) of the wave packet
     * @param u Vector variant of the matrix U
     */
    void set_U(const double x_c, const double y_c, const double sgm_x, const double sgm_y, const double p_x, const double p_y, arma::cx_vec &u_given);

    /**
     * @brief Evolve the system by a spatial step h and a time step dt
     * @param h Spatial step
     * @param dt Time step
     */
    void evolve(const double h, const double dt);
};

#endif //__schrodinger_hpp__