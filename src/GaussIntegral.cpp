//
// Created by Jakub Nowak on 21.10.2023.
//

#include "../include/GaussIntegral.h"

GaussIntegral::GaussIntegral(const int n){
    N = n;
    A = new double[N];
    x = new double[N];
    switch (N) {
        case 2:
            x[0] = -sqrt(1./3.);
            x[1] = sqrt(1./3.);
            A[0] = A[1] = 1.0;
            break;
        case 3:
            x[0] = -sqrt(3./5.);
            x[1] = 0.0;
            x[2] = sqrt(3./5.);
            A[0] = A[2] = 5.0 / 9.0;
            A[1] = 8.0 / 9.0;
            break;
        case 4:
            x[0] = -sqrt( 3./7. + 2./7.*sqrt(6./5.));
            x[1] = -sqrt( 3./7. - 2./7.*sqrt(6./5.));
            x[2] = sqrt( 3./7. - 2./7.*sqrt(6./5.));
            x[3] = sqrt( 3./7. + 2./7.*sqrt(6./5.));
            A[0] = A[3] = (18. - sqrt(30.)) / 36.;
            A[1] = A[2] = (18. + sqrt(30.)) / 36.;
            break;
        case 5:
            x[0] = 1./3. * sqrt(5. + 2. * sqrt(10./7.));
            x[1] = - 1./3. * sqrt(5. - 2. * sqrt(10./7.));
            x[2] = 0.0;
            x[3] = 1./3. * sqrt(5. - 2. * sqrt(10./7.));
            x[4] = - 1./3. * sqrt(5. + 2. * sqrt(10./7.));
            A[0] = A[4] = (322. - 13.*sqrt(70.)) / 900.;
            A[1] = A[3] = (322. + 13.*sqrt(70.)) / 900.;
            A[2] = 128. / 225.;
            break;
        default:
            cerr<<"\nIlosc punktow jest nieodpowiednia"<<endl;
            exit(-1);
    }
}
GaussIntegral::~GaussIntegral(){
    delete[] A;
    delete[] x;
}

double GaussIntegral::Gauss_1D(double(*f)(double)){
    result = 0.0;
    for(int i=0; i<this->N; i++){
        result += A[i] * f(x[i]);
    }
    return result;
}

double GaussIntegral::getXpc_dN(double(*f)(double), int i){
    return f(x[i]);
}

double GaussIntegral::getXpc_N(double(*f)(double, double), int i, int j){
    return f(x[i], x[j]);
}

double GaussIntegral::Gauss_2D(double(*f)(double, double)){
    result = 0.0;
    for(int i=0; i < this->N; i++){
        for (int j = 0; j < this->N; j++)
            result += A[i] * A[j] * f(x[i], x[j]);
    }
    return result;
}

double GaussIntegral::getA(const int i){
    return A[i];
}
double GaussIntegral::getX(const int i){
    return x[i];
}