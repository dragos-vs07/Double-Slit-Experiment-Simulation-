#pragma once
#include "simulation.h"
#include <cmath>
#include <cstdlib>

const int eta = 6, d = 10, cap_power = 2;

Simulation::Simulation() {}


complex_number gaussian_wave_packet_formula(double x,double y,const Params &_p)
{
    x *= _p.delta_x;
    y *= _p.delta_y;

    double a = - double(x-_p.x_zero)*double(x-_p.x_zero) / (4*pow(_p.sigma_x,2)) - double(y-_p.y_zero)*double(y-_p.y_zero) / (4*pow(_p.sigma_y,2)); // real number
    double b = _p.kx * double(x) + _p.ky * double(y); // real number

    double epow = std::exp(a);

    return complex_number( epow * cos(b), epow * sin(b) );  // complex number
}

void Simulation::normalise()
{
    double S = 0.0;

    for(int i=1;i<=NMAX-5;i++)
        for(int j=1;j<=MMAX-5;j++)
        {
            double abs_psi = psi_[i][j].absolute_value();
            S += abs_psi * abs_psi;
        }

    double A = 1.0 / sqrt(S);

    for(int i=1;i<=NMAX-5;i++)
        for(int j=1;j<=MMAX-5;j++)
            psi_[i][j] = psi_[i][j] * A;
}

double fx(double x)
{
    double xb = 10, xr = 140;
    if(x < xb)
        return pow((xb-x)/d,cap_power);
    else if(x > xr)
        return pow((x-xr)/30,cap_power);
    else
        return 0;

}

double fy(double y)
{
    double yb = 10, yr = 190;
    if(y < yb)
        return pow((yb-y)/d,cap_power);
    else if(y > yr)
        return pow((y-yr)/d,cap_power);
    else
        return 0;
}

void Simulation::initialise_V()
{
    int c1 = p_.slit_center - p_.slit_sep / 2, c2 = p_.slit_center + p_.slit_sep / 2;
    for (int i = 2; i <= 199; i++) {
        bool slit = abs(i - c1) <= p_.slit_w/2 || abs(i - c2) <= p_.slit_w/2;
        if (!slit)
            for (int j = 70; j < 74; j++) V_[i][j] = complex_number(20.0,0);
    }

    // x boundary ( left right )
    for(int i=1;i<=NMAX-5;i++)
    {
        for(int j=1;j<=9;j++)
            V_[i][j] = complex_number( V_[i][j].a, -eta*( fx(j) + fy(i) ) );

        for(int j=141;j<=MMAX-5;j++)
            V_[i][j] = complex_number( V_[i][j].a, -eta*( fx(j) + fy(i) ) );
    }

    // y boundary ( top bottom )
    for(int j=1;j<=MMAX-5;j++)
    {
        for(int i=1;i<=9;i++)
            V_[i][j] = complex_number( V_[i][j].a, -eta*( fx(j) + fy(i) ) );

        for(int i=191;i<=NMAX-5;i++)
            V_[i][j] = complex_number( V_[i][j].a, -eta*( fx(j) + fy(i) ) );
    }

}

void Simulation::reset(const Params &p)
{

    p_ = p;
    psi_.assign(NMAX, std::vector<complex_number>(MMAX));
    tmp_.assign(NMAX, std::vector<complex_number>(MMAX));
    K1_.assign(NMAX, std::vector<complex_number>(MMAX));
    K2_.assign(NMAX, std::vector<complex_number>(MMAX));
    K3_.assign(NMAX, std::vector<complex_number>(MMAX));
    K4_.assign(NMAX, std::vector<complex_number>(MMAX));
    V_.assign(NMAX, std::vector<complex_number>(MMAX));

    for(int i=1;i<=ny();i++)
        for(int j=1;j<=nx();j++)
            psi_[i][j] = gaussian_wave_packet_formula(j,i,p_);

    normalise();
    initialise_V();
    updateDensity();
}

void Simulation::step(int n)
{
    for(int ts=1;ts<=n;ts++)
    {
        deriv(psi_, K1_);

        for(int i=1;i<=NMAX-5;i++)
            for(int j=1;j<=MMAX-5;j++)
                tmp_[i][j] = psi_[i][j] + (p_.delta_t/2) * K1_[i][j];

        deriv(tmp_, K2_);

        for(int i=1;i<=NMAX-5;i++)
            for(int j=1;j<=MMAX-5;j++)
                tmp_[i][j] = psi_[i][j] + (p_.delta_t/2) * K2_[i][j];

        deriv(tmp_, K3_);

        for(int i=1;i<=NMAX-5;i++)
            for(int j=1;j<=MMAX-5;j++)
                tmp_[i][j] = psi_[i][j] + p_.delta_t * K3_[i][j];

        deriv(tmp_, K4_);

        for(int i=1;i<=NMAX-5;i++)
            for(int j=1;j<=MMAX-5;j++)
                psi_[i][j] = psi_[i][j] + p_.delta_t * ( K1_[i][j] + 2*K2_[i][j] + 2*K3_[i][j] + K4_[i][j] ) / 6;
    }
    updateDensity();
}
