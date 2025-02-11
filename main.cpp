// ------------------------------------------------------------------------
// | Metody Elementow Skonczonych, Jakub Nowak, semestr V, program glowny |
// ------------------------------------------------------------------------


#include "include/grid.h"
#include "include/GaussIntegral.h"
#include "include/UniElement.h"
#include "include/SoE.h"

const unsigned int gaussPoints = 2;     //punkty całkowania

//string path = "../input_files/Test1_4_4.txt";
//string path = "../input_files/Test2_4_4_MixGrid.txt";
string path = "../input_files/Test3_31_31_kwadrat.txt";
//string path = "../input_files/Test4_31_31_trapez.txt";
//string path = "../input_files/Test.txt";

using namespace std;

int main() {

    //------------------------------------------ DANE ------------------------------------------
    //do obliczen czasu
    time_t start, end, tResult;
    time(&start);           //start pomiaru czasu

    //plik do odczytu
    ifstream inputFile(path.c_str());
    if (!inputFile.good()) {
        cerr << "Blad otwarcia pliku" << endl;
    }

    //siatka
    Grid grid;
    GlobalData globalData;

    //wczytywanie z pliku (wszystkiego)
    readFromFile(inputFile, globalData, grid);
    inputFile.close();

    //wyswietlanie danych ogolnych
    //globalData.printGlobalData();

    //element uniwersalny dla okreslonej ilosci pkt calkownia
    UniElement uElement(gaussPoints);
    //uElement.printUniElements();

    //------------------------------------------ OBLICZENIA ------------------------------------------
    //obliczanie macierzy H i Hbc - warunek brzegowy
    try{
        for (int i = 0; i < globalData.Elements_number; ++i){
            uElement.calculateH_C(grid.elements[i], globalData);
            uElement.calculateHbc_P(grid.elements[i], globalData);
            grid.elements[i].calculateGlobalH(globalData);
        }

        //obliczanie ukladu rownan
        SoE systemOfEq(globalData);

        //agregacja
        for (int i = 0; i < globalData.Elements_number; ++i){
            systemOfEq.aggregation(grid.elements[i]);
        }

        //uklad rownan
        try{
            //obiczenie i zapis do pliku!
            systemOfEq.calcuateResult(globalData, grid);
        }
        catch(string e) {
            cerr<<"Uklad rownan --> "<<e<<endl;
        }

    }catch (string e){
        cerr<<"Obliczenia: "<<e<<endl;
    }

    time(&end);
    tResult = end - start;
    cout<<endl<<"Czas obliczen: "<<tResult<<" s"<<endl;
    //wyswietlenie siatki
    //grid.printGrid();

    return 0;
}

