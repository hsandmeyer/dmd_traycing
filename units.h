#include <cmath>
#ifndef UNITS_H
#define UNITS_H
constexpr long double operator"" _cm(long double x){
    return x / 100.;
}

constexpr long double operator"" _mm(long double x){
    return x / 1000.;
}

constexpr long double operator"" _nm(long double x){
    return x * 1.e-9;
}

constexpr long double operator"" _um(long double x){
    return x * 1.e-6;
}

constexpr long double operator"" _deg(long double x){
    return x / 180. * M_PI;
}

constexpr double operator"" _deg(unsigned long long x){
    return x / 180. * M_PI;
}
#endif
