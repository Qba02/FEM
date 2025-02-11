// Created by Jakub Nowak on 21.10.2023.
// --------------------------------------------------
// | Klasa do obliczenia całki metoda Gauusa-Crouta |
// --------------------------------------------------

#include "grid.h"

#ifndef PROJEKT_MES_GAUSSINTEGRAL_H
#define PROJEKT_MES_GAUSSINTEGRAL_H

using namespace std;

class GaussIntegral{
    double result = 0.0;
    int N = 0;

protected:
    double *x;
    double *A;

public:
    explicit GaussIntegral(const int n);
    GaussIntegral()= default;
    ~GaussIntegral();

    double Gauss_1D(double(*f)(double));
    double Gauss_2D(double(*f)(double, double));

    double getXpc_dN(double(*f)(double), int i);
    double getXpc_N(double(*f)(double, double), int i, int j);
    double getA(const int i);
    double getX(const int i);
};


#endif //PROJEKT_MES_GAUSSINTEGRAL_H
