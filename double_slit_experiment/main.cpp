#include <bits/stdc++.h>
using namespace std;

ofstream fout("output.txt");

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

double V[205][125]; // V[i][j] = potential energy at cell (i,j)

class complex_number{   // defining complex numbers

    private:

    double a,b;

    public:

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

complex_number M[205][125][105];  // M[i][j][n] = wave function at cell (i,j) and time stamp n
complex_number L[205][125]; // L[i][j] = Laplacian at cell (i,j)
complex_number complex_i = complex_number(0,1); // just a constant

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
void initialise_V()
{
    int slit_w = 5, c1 = 80, c2 = 120;
    for (int i = 1; i <= 200; i++) {
        bool slit = abs(i - c1) <= slit_w/2 || abs(i - c2) <= slit_w/2;
        if (!slit)
            for (int j = 70; j < 74; j++) V[i][j] = 20.0;
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
void calculate_laplacian(int n)
{
    for(int i=1;i<=200;i++)
        for(int j=1;j<=120;j++)
        {
            complex_number lap_x = ( M[i][j-1][n] + M[i][j+1][n] - 2 * M[i][j][n] ) / ( delta_x * delta_x );
            complex_number lap_y = ( M[i-1][j][n] + M[i+1][j][n] - 2 * M[i][j][n] ) / ( delta_y * delta_y ) ;

            L[i][j] = lap_x + lap_y;
        }
}
void evolve(int n)
{

    for(int ts=1;ts<=n;ts++)
    {
        calculate_laplacian(ts-1);
        for(int i=1;i<=200;i++)
            for(int j=1;j<=120;j++)
            {
                complex_number term1 = ( complex_i * h ) / ( 2 * e_mass ); // complex number
                complex_number term2 = complex_i / h; // complex number

                M[i][j][ts] = M[i][j][ts-1] + delta_t * ( term1 * L[i][j] - term2 * V[i][j] * M[i][j][ts-1] );
            }
    }
}
void dump_probability_results(int n)
{
    // dumps in a text file the final probability matrix
    // results should be taken to python to plot
    for(int i=1;i<=200;i++)
    {
        for(int j=1;j<=120;j++)
        {
            double abs_psi = M[i][j][n].absolute_value();
            fout << abs_psi * abs_psi <<" ";
        }
        fout << '\n';
    }
}

int main()
{
    int n = 100;
    first_timestamp();
    initialise_V(); // potential energy doesn't take part in defining the starting phase of the packet
    calculate_A();
    normalise();
    evolve(n);    // takes in number of steps
    dump_probability_results(n);

    return 0;
}
