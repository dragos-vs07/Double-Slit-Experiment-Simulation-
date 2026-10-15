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

complex_number M[205][125][110];  // M[i][j][n] = wave function at cell (i,j) and time stamp n
complex_number cur[205][125], tmp[205][125], K1[205][125], K2[205][125], K3[205][125], K4[205][125]; 
complex_number complex_i = complex_number(0,1); // just a constant
complex_number V[205][125]; // V[i][j] = CAP

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
        for(int j=1;j<=120;j++)
        M[i][j][0] = gaussian_wave_packet_formula(j,i);
}
double n = 2, eta = 6, d = 10;
double fx(double x)
{
    double xb = 10, xr = 110;
    if(x < xb)
        return pow((xb-x)/d,n);
    else if(x > xr)
        return pow((x-xr)/d,n);
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

        for(int j=111;j<=120;j++)
        V[i][j] = complex_number( V[i][j].a, -eta*( fx(j) + fy(i) ) );
    }

    // y boundary ( top bottom )
    for(int j=1;j<=120;j++)
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
        for(int j=1;j<=120;j++)
        {
            double abs_psi = M[i][j][0].absolute_value();
            S += abs_psi * abs_psi;
        }

    A = 1.0 / sqrt(S);
}
void normalise()
{
    for(int i=1;i<=200;i++)
        for(int j=1;j<=120;j++)
            M[i][j][0] = M[i][j][0] * A;
}
void deriv(complex_number inM[205][125], complex_number outM[205][125])
{
    complex_number term1 = (complex_i * h) / (2 * e_mass);
    complex_number term2 = complex_i / h;

    for(int i=1;i<=200;i++)
        for(int j=1;j<=120;j++)
        {
            complex_number lap_x = ( inM[i][j-1] + inM[i][j+1] - 2 * inM[i][j] ) / ( delta_x * delta_x );
            complex_number lap_y = ( inM[i-1][j] + inM[i+1][j] - 2 * inM[i][j] ) / ( delta_y * delta_y ) ;

            outM[i][j] = term1 * (lap_x + lap_y) - term2 * V[i][j] * inM[i][j];
        }
}
void evolve(int n)
{
    for(int i=1;i<=200;i++)
        for(int j=1;j<=120;j++)
            cur[i][j] = M[i][j][0];

    for(int ts=1;ts<=n;ts++)
    {
        deriv(cur, K1);

        for(int i=1;i<=200;i++)
            for(int j=1;j<=120;j++)
                tmp[i][j] = cur[i][j] + (delta_t/2) * K1[i][j];
        deriv(tmp, K2);

        for(int i=1;i<=200;i++)
            for(int j=1;j<=120;j++)
                tmp[i][j] = cur[i][j] + (delta_t/2) * K2[i][j];
        deriv(tmp, K3);

        for(int i=1;i<=200;i++)
            for(int j=1;j<=120;j++)
                tmp[i][j] = cur[i][j] + delta_t * K3[i][j];
        deriv(tmp, K4);

        for(int i=1;i<=200;i++)
            for(int j=1;j<=120;j++)
                cur[i][j] = cur[i][j]
                    + (delta_t/6) * (K1[i][j] + 2*K2[i][j] + 2*K3[i][j] + K4[i][j]);

        if(ts % 10 == 0)
        {
            double S = 0;
            for(int i=1;i<=200;i++)
                for(int j=1;j<=120;j++)
                {
                    M[i][j][ts/10] = cur[i][j];
                    double a = cur[i][j].absolute_value();
                    S += a * a;
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
            for(int j=1;j<=120;j++)
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
    int n = 1050;
    first_timestamp();
    initialise_V(); // potential energy doesn't take part in defining the starting phase of the packet
    calculate_A();
    normalise();
    evolve(n);    // takes in number of steps
    dump_probability_results(n);

    return 0;
}
