// Created by Jakub Nowak on 29.11.2023.
// ------------------------------------------
// | Klasa do rozwiazaywania ukladow rownan |
// ------------------------------------------
#include "grid.h"

#ifndef PROJEKT_MES_SOE_H
#define PROJEKT_MES_SOE_H



class SoE {     //system of equations - układ rownan
private:
    matrix globalH;
    matrix globalP;
    matrix globalC;
    matrix X;             //rozwiazanie ukladu rownan
    vector<double> resultMax;     // wektor wynikow max
    vector<double> resultMin;     // wektor wynikow min

public:
    static double eps;     //dokladnosc
    static int nodesN;    //ilosc wezlow - wymiar macierzy - liczba rownan do rozwiazania

    explicit SoE(GlobalData& gd);
    void aggregation(const Element& e);
    void Gauss_Crout();
    void printAndAgregate();
    void calcuateResult(const GlobalData& gd, const Grid& grid);
    void test(double *tMax, double *tMin);
};


#endif //PROJEKT_MES_SOE_H
