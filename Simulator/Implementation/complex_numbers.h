#ifndef COMPLEX_NUMBERS_H
#define COMPLEX_NUMBERS_H

#pragma once

#include <cmath>

class complex_number
{
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

#endif // COMPLEX_NUMBERS_H
