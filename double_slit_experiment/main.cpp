#include <bits/stdc++.h>
using namespace std;

// all over the code when we use cell coordinates in formulas we pass it as x = j and y = i because the x axis is the columns and y axis the lines

double e = 2.718281828459045; // euler's constant
double h = 1; // reduced plank
double e_mass = 1; // mass of electron
double delta_x = 0.2, delta_y = 0.2;    // real distance between cells
double delta_t = 0.005; // time per timestamp
double kx = 3, ky = 0;  // gaussian packet momentum
double x_zero = 6.0,y_zero = 20.0;  // origin of the packet ( position of electron gun )
double sigma_x = 2, sigma_y = 4; // dimensions of the packet
double A; // normalisation constant

class complex_number{   // defining complex numbers

    public:

    double a,b;
    // constructors

    complex_number(double x, double y)
    {
        a = x;
        b = y;
    }

    complex_number()
    {
        a = 0;
        b = 0;
    }

    // addition
    friend complex_number operator +(complex_number c1, complex_number c2)
    {
        return complex_number(c1.a + c2.a, c1.b + c2.b);
    }
    friend complex_number operator +(double c1, complex_number c2)
    {
        return complex_number(c1 + c2.a, c2.b);
    }
    friend complex_number operator +(complex_number c1, double c2)
    {
        return complex_number( c1.a + c2, c1.b );
    }

    // subtraction
    friend complex_number operator -(complex_number c1, complex_number c2)
    {
        return complex_number( c1.a - c2.a, c1.b - c2.b );
    }
    friend complex_number operator -(double c1, complex_number c2)
    {
        return complex_number( c1 - c2.a, -c2.b );
    }
    friend complex_number operator -(complex_number c1, double c2)
    {
        return complex_number( c1.a - c2, c1.b );
    }

    // multiplication
    friend complex_number operator *(complex_number c1, double c2)
    {
        return complex_number( c1.a * c2, c1.b * c2 );
    }
    friend complex_number operator *(double c1, complex_number c2)
    {
        return complex_number( c2.a * c1, c2.b * c1 );
    }
    friend complex_number operator *(complex_number c1,complex_number c2)
    {
        return complex_number( c1.a * c2.a - c1.b * c2.b, c1.a * c2.b + c1.b * c2.a );
    }

    // division by a real
    friend complex_number operator /(complex_number c1, double c2)
    {
        return complex_number( c1.a / c2, c1.b / c2 );
    }

    // absolute value
    double absolute_value()
    {
        return sqrt(a * a + b * b);
    }
};

const int NMAX = 205, YMAX = 175;

complex_number M[NMAX][YMAX][130];  // M[i][j][n] = wave function at cell (i,j) in the n-th snapshot
complex_number psi[NMAX][YMAX], tmp[NMAX][YMAX], K1[NMAX][YMAX], K2[NMAX][YMAX], K3[NMAX][YMAX], K4[NMAX][YMAX]; 
complex_number complex_i = complex_number(0,1); // just a constant
complex_number V[NMAX][YMAX]; // V[i][j] = potential energy ( cap near borders )

complex_number gaussian_wave_packet_formula(double x,double y)
{
    x *= delta_x;
    y *= delta_y;

    double a = - double(x-x_zero)*double(x-x_zero) / (4*pow(sigma_x,2)) - double(y-y_zero)*double(y-y_zero) / (4*pow(sigma_y,2)); // real number
    double b = kx * double(x) + ky * double(y); // real number

    double epow = pow(e,a);

    return complex_number( epow * cos(b), epow * sin(b) );  // complex number
}

void first_timestamp()
{
    for(int i=1;i<=200;i++)
        for(int j=1;j<=170;j++)
        M[i][j][0] = gaussian_wave_packet_formula(j,i);
}
double n = 2, eta = 6, d = 10;
double fx(double x)
{
    double xb = 10, xr = 140;
    if(x < xb)
        return pow((xb-x)/d,n);
    else if(x > xr)
        return pow((x-xr)/30,n);
    else
        return 0;

}
double fy(double y)
{
    double yb = 10, yr = 190;
    if(y < yb)
        return pow((yb-y)/d,n);
    else if(y > yr)
        return pow((y-yr)/d,n);
    else
        return 0;
}
void initialise_V()
{
    int slit_w = 5, c1 = 80, c2 = 120;
    for (int i = 2; i <= 199; i++) {
        bool slit = abs(i - c1) <= slit_w/2 || abs(i - c2) <= slit_w/2;
        if (!slit)
            for (int j = 70; j < 74; j++) V[i][j] = complex_number(20.0,0);
    }

    // x boundary ( left right )
    for(int i=1;i<=200;i++)
    {
        for(int j=1;j<=9;j++)
        V[i][j] = complex_number( V[i][j].a, -eta*( fx(j) + fy(i) ) );

        for(int j=141;j<=170;j++)
        V[i][j] = complex_number( V[i][j].a, -eta*( fx(j) + fy(i) ) );
    }

    // y boundary ( top bottom )
    for(int j=1;j<=170;j++)
    {
        for(int i=1;i<=9;i++)
        V[i][j] = complex_number( V[i][j].a, -eta*( fx(j) + fy(i) ) );

        for(int i=191;i<=200;i++)
        V[i][j] = complex_number( V[i][j].a, -eta*( fx(j) + fy(i) ) );
    }
        
}
void calculate_A() 
{
    double S = 0.0;

    for(int i=1;i<=200;i++)
        for(int j=1;j<=170;j++)
        {
            double abs_psi = M[i][j][0].absolute_value();
            S += abs_psi * abs_psi;
        }

    A = 1.0 / sqrt(S);
}
void normalise()
{
    for(int i=1;i<=200;i++)
        for(int j=1;j<=170;j++)
            M[i][j][0] = M[i][j][0] * A;
}
void deriv(complex_number psi[NMAX][YMAX], complex_number psi_deriv[NMAX][YMAX])
{
    // calculate derivative of time of psi by Schrodinger formula for a known psi ( all cells )

    complex_number term1 = (complex_i * h) / (2 * e_mass);
    complex_number term2 = complex_i / h;

    for(int i=1;i<=200;i++)
        for(int j=1;j<=170;j++)
        {
            complex_number lap_x = ( psi[i][j-1] + psi[i][j+1] - 2 * psi[i][j] ) / ( delta_x * delta_x );
            complex_number lap_y = ( psi[i-1][j] + psi[i+1][j] - 2 * psi[i][j] ) / ( delta_y * delta_y ) ;
            // lap_x + lap_y = Laplacian

            psi_deriv[i][j] = term1 * (lap_x + lap_y) - term2 * V[i][j] * psi[i][j];
        }
}
void evolve(int n)
{
    // psi at ts = 0 holds the wave function for every cell calculated with the gaussian package
    // which is also always the first snapshot

    for(int i=1;i<=200;i++)
        for(int j=1;j<=170;j++)
            psi[i][j] = M[i][j][0];

    for(int ts=1;ts<=n;ts++)
    {
        // at each time step we first find psi(t) differentiated with respect to t for every cell

        deriv(psi, K1); // K1, simple Euler, finds slope of psi at time point ts for every cell
        
        // now in temp we hold all the points from psi moved by delta_t / 2 time units forward

        for(int i=1;i<=200;i++)
            for(int j=1;j<=170;j++)
                tmp[i][j] = psi[i][j] + (delta_t/2) * K1[i][j]; // slope at ( ts + delta_t ) / 2
                                                                // by direction of calculated slope K1

        // so we calculate the slope of psi at time point delta_t / 2 along the direction K1 slopes indicate
        deriv(tmp, K2); // K2, finds the slope of all new calculated points after the delta_t / 2 time movement

        // repeat the same process but instead of moving on the direction of K1 ( so slope at delta_t + ts)
        // we move along the direction of K2 ( slope at delta_t/2 + ts)

        for(int i=1;i<=200;i++)
            for(int j=1;j<=170;j++)
                tmp[i][j] = psi[i][j] + (delta_t/2) * K2[i][j]; // slope at ts + delta_t / 2
                                                                // by direction of calculated slope K2

        deriv(tmp, K3); // more accurate slope for all points of psi at time moment ts + delta_t/2

        // now temp holds all the points of psi at the next time step ( so ts + delta_t )

        for(int i=1;i<=200;i++)
            for(int j=1;j<=170;j++)
                tmp[i][j] = psi[i][j] + delta_t * K3[i][j]; // slope at ts + delta_t


        deriv(tmp, K4); // and so we calculate the slope of all the cells of psi at the next time step ( ts + delta_t) 

        for(int i=1;i<=200;i++)
            for(int j=1;j<=170;j++)
                psi[i][j] = psi[i][j] + delta_t * ( K1[i][j] + 2*K2[i][j] + 2*K3[i][j] + K4[i][j] ) / 6; 

            // in the end, RK4 approximation for psi differentiated with respect to t is the avrage of all 4 slopes calculated:
            // the one at time moment ts , two at time moment ts + delta_t/2 ( most important so multiplied by 2 ) and one at the next time step ts + delta_t
            // thus psi at the next timestep is line 248

        if(ts % 10 == 0) // save snapshot at every 10 time steps
        {
            double S = 0;
            for(int i=1;i<=200;i++)
                for(int j=1;j<=170;j++)
                {
                    M[i][j][ts/10] = psi[i][j]; // store snapshot
                    double a = psi[i][j].absolute_value();
                    S += a * a; // sum to check probability correctness
                }
            //cout << ts << " " << S << "\n";   // should stay ~1
        }
    }
}
void dump_probability_results(int n)
{
    // dumps in more text files "snapshots" of the probability matrix
    // results should be taken to python to plot

    for(int s=0;s<=n;s+=10) // incrementation should be modified if n is also highly modified 
    {
        ofstream fout( "C:/Users/Dragos/Desktop/Double Slit Experiment Simulation/double_slit_experiment/snapshots/output_" + to_string(s) + ".txt" );
        for(int i=1;i<=200;i++)
        {
            for(int j=1;j<=170;j++)
            {
                double abs_psi = M[i][j][s/10].absolute_value();
                fout << abs_psi * abs_psi <<" ";
            }
            fout << '\n';
        }
    }
}

int main()
{
    int n = 1200; // time steps
    first_timestamp();
    initialise_V(); // potential energy doesn't take part in defining the starting phase of the packet
    calculate_A();
    normalise();
    evolve(n);    // takes in number of steps
    dump_probability_results(n);

    return 0;
}
