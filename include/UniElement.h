// Created by Jakub Nowak on 08.11.2023.
// -----------------------------------------------------------------
// | Element Uniwersalny - obliczanie macierzy H, C oraz wektora P |
// -----------------------------------------------------------------


#include "grid.h"
#include "GaussIntegral.h"
#include "matrix.h"
#ifndef PROJEKT_MES_UNIELEMENT_H
#define PROJEKT_MES_UNIELEMENT_H

class Jacobian {
public:
    matrix J;
    matrix Jinv;    //macierz odwrotna
    double detJ = 0.;

    Jacobian(int N, int M): J(N, M), Jinv(N, M){}

};


class UniElement{
private:
    matrix ksi_dN;   //mcierze rozmiaru dN x node, gdzie node to liczba wezlow a dN to zawsze 4 (w 2D) dla macierzy H!! - pochodne
    matrix eta_dN;
    matrix ksiEta_N;           //macierz rozmiaru dN x node dla macierzy C!! - funkcje kształtu
    static const int dim = 2;    //wymiar - 2D
    static const int dN = 4;    //liczba wezlow elementu = liczba pochodnych w tym wymiarze, czyli liczba kolumn dksi/deta
    int Npc;                    //punkty calkowania
    int node;                   //ilosc punktow calkowania^2 (w praktyce liczba wierszy dla wielu macierzy)

    GaussIntegral integral;     //do liczenia calek
    vector <matrix> surface;    //boki kwadratu
public:
    explicit UniElement(int npc);
    void calculateH_C(Element& e, const GlobalData& gd);
    void calculateHbc_P(Element& e, const GlobalData& gd);
    void printUniElements() const;
    void calculateKsiEta();
    void calculateSurfaceN();       //czyli tablice dla kazdego boku
};

//pochodne funkcji kształtu --> n to eta_dN, e to ksi_dN
static double d1de(double n){return (-1./4.) * (1. - n);}
static double d2de(double n){return (1./4.) * (1. - n);}
static double d3de(double n){return (1./4.) * (1. + n);}
static double d4de(double n){return (-1./4.) * (1. + n);}
static double d1dn(double e){return (-1./4.) * (1. - e);}
static double d2dn(double e){return (-1./4.) * (1. + e);}
static double d3dn(double e){return (1./4.) * (1. + e);}
static double d4dn(double e){return (1./4.) * (1. - e);}

//funkcje kształtu
static double N1(double e, double n){return (1./4.) * (1. - e) * (1 - n);}
static double N2(double e, double n){return (1./4.) * (1. + e) * (1 - n);}
static double N3(double e, double n){return (1./4.) * (1. + e) * (1 + n);}
static double N4(double e, double n){return (1./4.) * (1. - e) * (1 + n);}

#endif //PROJEKT_MES_UNIELEMENT_H
