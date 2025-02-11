//
// Created by Jakub Nowak on 29.11.2023.
//

#include "../include/SoE.h"
#include <iomanip>

// -----------------------------------------------------------
// | UKLAD ROWNAN: ([H] * [C]/dT)*{t1} - ([C]/dT)*{t0} + {P} |
// -----------------------------------------------------------
SoE::SoE(GlobalData& gd): globalH(gd.Nodes_number, gd.Nodes_number), globalC(gd.Nodes_number, gd.Nodes_number),
                          globalP(gd.Nodes_number, 1), X(gd.Nodes_number, 1){
    X = matrix(gd.Nodes_number, 1, gd.InitialTemp);                    // wektor pionowy t0 - startowy
    nodesN = gd.Nodes_number;
}

double SoE::eps = 1e-12;
int SoE::nodesN;

void SoE::Gauss_Crout()
{
    matrix AB(nodesN, nodesN + 1);   //macierz do obliczen

    for(int i = 0; i < nodesN ; i++)
        for(int k = 0; k <= nodesN ; k++)
            if(k == nodesN)
                AB(i, k) = X(i, 0);
            else
                AB(i, k) = globalH(i, k);

    //pomocnicza tablica zamiany kolumn
    int *W = new int[nodesN];
    for(int i=0; i < nodesN; i++)
        W[i] = i;

    //iteratory petli
    int    i, j, k;
    double m, s;        //zmienne pomocnicze

    for( i = 0; i <= nodesN; i++ )
        W[i] = i;

    // eliminacja współczynników
    for( i = 0; i < nodesN - 1; i++ )
    {
        k = i;
        for( j = i + 1; j < nodesN; j++ )
            if( abs(AB(i, W[k])) < abs(AB(i, W[j])))
                k = j;
        swap (W[k], W[i]);


        for( j = i + 1; j < nodesN ; j++ )
        {
            if( abs(AB(i, W[i])) < eps ) throw string("Gauss_Crout:\nblad w rozwiazywaniu ukladu rownan (wspolczynniki)");
            m = - AB(j, W[i]) / AB(i, W[i]) ;
            for( k = i + 1; k <= nodesN; k++ )
                AB(j, W[k]) += m * AB(i,W[k]);
        }
    }

    // wyliczanie niewiadomych
    for( i = nodesN - 1; i >= 0; i--)
    {
        if( abs(AB(i, W[i])) < eps ) throw string("Gauss_Crout:\nblad w rozwiazywaniu ukladu rownan (niewiadome)");
        s = AB(i, nodesN);
        for( j = nodesN - 1; j >= i + 1; j--)
            s -= AB(i, W[j]) * X(W[j], 0);
        X(W[i], 0) = s / AB(i, W[i]);
    }

}

void SoE::aggregation(const Element &e) {

    int* n = get_size(e.H); //n[0] - liczba wierszy, n[1] - kolumn

    for (int i = 0; i < n[0]; ++i) {
        for (int j = 0; j < n[1]; ++j) {
            globalH(e.nodes[i].node_id-1, e.nodes[j].node_id-1) += e.H(i, j);
            globalC(e.nodes[i].node_id-1, e.nodes[j].node_id-1) += e.C(i, j);
        }
        globalP(e.nodes[i].node_id-1, 0) += e.P(i, 0);
    }
}

void SoE::calcuateResult(const GlobalData& gd, const Grid& grid) {
    int symStep = (int)gd.SimulationStepTime;
    globalC = globalC/symStep;                        // macierz [C]/dT              //T=tał(czas)
    globalH = globalH + globalC;                      // ([H] * [C]/dT) * {t1} - układ rownan


    for (int i = symStep; i <= gd.SimulationTime; i+=symStep) {
        ostringstream filename;
        filename << "../output_files/Test3/foo" << i << ".vtk";

        //plik do zapisu
        ofstream outputFile(filename.str());
        if (!outputFile.good()) {
            cerr <<"Blad otwarcia pliku" << endl;
            exit(0);
        }

        X = globalC * X + globalP;                   // ([C]/dT)*{t0} + {P} - wyrazy wolne
        this->Gauss_Crout();
        writeToFile(outputFile, gd, grid, X);

        //wydruk calego wyniku lub samych wartosci max/min
        cout<<endl<<"t: "<<i<<"s "<<endl;
        this->printAndAgregate();
        outputFile.close();
    }

}

void SoE::printAndAgregate() {
//    int * n = get_size(globalC);
//    cout<<"Globalna macierz C"<<endl;
//    for (int i = 0; i < n[0]; ++i) {
//        for (int j = 0; j < n[1]; ++j) {
//            cout<<setprecision(3)<<setw(5)<<globalC(i, j)<<"|";
//        }
//        cout<<endl;
//    }
//    cout<<endl<<"Globalny wektor P"<<endl;
//    for (int i = 0; i < n[0]; ++i) {
//        cout<<setprecision(7)<<globalP(i, 0)<<" | ";
//    }
//    cout<<trans(X)<<endl;

    int* n = get_size(X);
    double max=X(0, 0), min=X(0,0);
    for (int i = 0; i < n[0]; ++i) {
        if(max<X(i,0)){
            max = X(i,0);
        }else if(min>X(i,0)){
            min = X(i,0);
        }
    }
    resultMax.push_back(max);
    resultMin.push_back(min);
}

void SoE::test(double *tMax, double *tMin){
    //test
    int* n = get_size(X);
    double epsilon = 1e-01;
    cout<<"-----------------"<<"TESTOWANIE Z DOKLADNOSCIA: "<<epsilon<<"-----------------"<<endl;
    cout<<setw(13)<<"PROGRAM "<<setw(10)<<"TEST"<<setw(22)<<"PROGRAM "<<setw(10)<<"TEST"<<endl;
    for (int i = 0; i < 10; ++i) {      // testowanie 10 pierwszych wartosci
        if(abs(resultMin[i] - tMin[i]) > epsilon || abs(resultMax[i] - tMax[i]) > epsilon ){
            cout<<resultMin[i]<<" - "<<tMin[i]<<" = "<<abs(resultMin[i] - tMin[i])<<endl;
            cout<<resultMax[i]<<" - "<<tMax[i]<<" = "<<abs(resultMax[i] - tMax[i])<<endl;
            cerr<<"Blad - wyniki nie sa zgodne"<<endl;
            exit(0);
        }else{
            cout<<setprecision(10)<<"Min: "<<resultMin[i]<<" | "<<tMin[i]<<"\tMax: "<<tMax[i]<<" | "<<resultMax[i]<<endl;
        }
    }
    cout<<endl<<"Test OK"<<endl;
}