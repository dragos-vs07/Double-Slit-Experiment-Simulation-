#ifndef SIMULATION_H
#define SIMULATION_H

#pragma once

#include <vector>
#include "complex_numbers.h"

struct Params{
    double delta_x = 0.2, delta_y = 0.2;    // real distance between cells
    double delta_t = 0.005; // time per timestamp
    double kx = 3, ky = 0;  // gaussian packet momentum
    double x_zero = 6.0,y_zero = 20.0;  // origin of the packet ( position of electron gun )
    double sigma_x = 2, sigma_y = 4; // dimensions of the packet

    int slit_w = 5;        // slit width (cells)
    int slit_center = 100; // midpoint between the two slits (y, in cells)
    int slit_sep = 40;     // distance between the slit centers (your c1=80, c2=120)
};

using Grid = std::vector<std::vector<complex_number>>;

class Simulation
{
public:
    Simulation();
    void reset(const Params &p);   // build starting state from p
    void normalise();
    void initialise_V();
    void step(int n);              // advance n RK4 steps

    const std::vector<double>& density() const
    {
        return density_;
    }

    int nx() const
    {
        return MMAX - 5; // 170 columns (x)
    }
    int ny() const
    {
        return NMAX - 5; // 200 rows (y)
    }

private:
    static constexpr int NMAX = 205, MMAX = 175;

    Params p_;                               // copy of the current parameters
    Grid psi_, tmp_, K1_, K2_, K3_, K4_, V_;

    complex_number complex_i = complex_number(0,1); // constants
    double e = 2.718281828459045; // euler's constant
    double h = 1; // reduced plank
    double e_mass = 1; // mass of electron

    std::vector<double> density_;

    void deriv(const Grid &psi, Grid &psi_deriv) // calculate derivative of time of psi by Schrodinger formula for a known psi ( all cells )
    {

        complex_number term1 = (complex_i * h) / (2 * e_mass);
        complex_number term2 = complex_i / h;

        for(int i=1;i<=NMAX-5;i++)
            for(int j=1;j<=MMAX-5;j++)
            {
                complex_number lap_x = ( psi[i][j-1] + psi[i][j+1] - 2 * psi[i][j] ) / ( p_.delta_x * p_.delta_x );
                complex_number lap_y = ( psi[i-1][j] + psi[i+1][j] - 2 * psi[i][j] ) / ( p_.delta_y * p_.delta_y ) ;
                // lap_x + lap_y = Laplacian

                psi_deriv[i][j] = term1 * (lap_x + lap_y) - term2 * V_[i][j] * psi[i][j];
            }
    }
    void updateDensity()   // fill density_ from psi_
    {
        int k = -1;
        density_.resize(nx()*ny());
        for(int i=1;i<=ny();i++)
            for(int j=1;j<=nx();j++)
            {
                double abs_psi = psi_[i][j].absolute_value();
                density_[++k] = abs_psi * abs_psi;
            }

    }

};

#endif // SIMULATION_H
