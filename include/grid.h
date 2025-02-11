// Created by Jakub Nowak on 18.10.2023.
// ------------------------------------------------------------------
// | Zmienne globalne, struktura siatki, wczytywanie/zapis do pliku |
// ------------------------------------------------------------------

#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <iomanip>
#include <vector>
#include "matrix.h"
#include <ctime>

#ifndef PROJEKT_MES_DATA_H
#define PROJEKT_MES_DATA_H

using namespace std;

//-------------------------------------GlobaData----------------------------------
class GlobalData{
public:
    double SimulationTime, SimulationStepTime, Conductivity,
            Alfa, Tot, InitialTemp, Density, SpecificHeat;
    int Nodes_number, Elements_number;
    GlobalData();
    void printGlobalData() const;
};

//-------------------------------------Node---------------------------------------
class Node{
public:
    double x,y;
    int node_id;
    short BC=0;       //flaga dla warunku brzegowego

    Node();
    Node(int ID, double wsp_x, double wsp_y);
    friend ostream& operator<<(ostream& output, const Node& obiekt);
};

//-------------------------------------Element-------------------------------------
class Element{
public:
    const static int node_number=4; //iosc wezlow w figurze
    int element_id;
    Node nodes[node_number];
    matrix H;       //macierz 4 x 4 -> potem dodajemy do niej Hbc
    matrix Hbc;     //macierze 4 x 4
    matrix P;       //macierze 4 x 1
    matrix C;       //macierze 4 x 4

    Element();
    Element(int ID, Node ur, Node ul, Node dl, Node dr);
    friend ostream& operator<<(ostream& output, const Element& obiekt);
    void calculateGlobalH(const GlobalData& gd);
};

//-------------------------------------Grid---------------------------------------
class Grid{
public:
    vector <Node> nodes;
    vector <Element> elements;

    void printGrid();
    friend void readFromFile(ifstream& file, GlobalData& gd, Grid& grid);
};

void writeToFile(ofstream& file, const GlobalData& gd, const Grid& grid, const matrix& X);     //zapis rezultatu

#endif //PROJEKT_MES_DATA_H
