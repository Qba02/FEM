//
// Created by Jakub Nowak on 08.11.2023.
//

#include "../include/UniElement.h"

UniElement::UniElement(int npc): integral(npc), ksi_dN((int)pow(npc, dim), dN),
                                 eta_dN((int)pow(npc, dim), dN), ksiEta_N((int)pow(npc, dim), dN){
    Npc = npc;
    node = (int)pow(npc, dim);    //tyle wierszy ile punktow calkowania
    surface.assign(dN, matrix(Npc, dN));
    calculateKsiEta();
    calculateSurfaceN();
}

void UniElement::calculateKsiEta() {
    int e = -1;         //e i n konieczne do Gaussa -> opowiadają indeksom wag/wartosci
    int n = -1;
    for (int j = 0; j < node; j++) {

        if(e == (Npc - 1)) e = 0;       // 0, 1, 0, 1
        else e++;

        if(j%Npc == 0) n++;           // 0, 0, 1, 1

        ksi_dN(j, 0) = integral.getXpc_dN(d1de, n);
        ksi_dN(j, 1) = integral.getXpc_dN(d2de, n);
        ksi_dN(j, 2) = integral.getXpc_dN(d3de, n);
        ksi_dN(j, 3) = integral.getXpc_dN(d4de, n);

        eta_dN(j, 0) = integral.getXpc_dN(d1dn, e);
        eta_dN(j, 1) = integral.getXpc_dN(d2dn, e);
        eta_dN(j, 2) = integral.getXpc_dN(d3dn, e);
        eta_dN(j, 3) = integral.getXpc_dN(d4dn, e);

        //cout<<integral.getX(e)<<" "<<integral.getX(n)<<endl;
        ksiEta_N(j, 0) = integral.getXpc_N(N1, e, n);
        ksiEta_N(j, 1) = integral.getXpc_N(N2, e, n);
        ksiEta_N(j, 2) = integral.getXpc_N(N3, e, n);
        ksiEta_N(j, 3) = integral.getXpc_N(N4, e, n);
    }
}

void UniElement::calculateSurfaceN() {

    double e = 0;         //e - ksi, n - eta
    double n = 0;

    for (int k = 0; k < dN; ++k) {

        for (int j = 0; j < Npc; ++j) {
            switch(k){
                case 0: e = integral.getX(j); n = -1.;
                    break;
                case 1: e = 1.; n = integral.getX(j);
                    break;
                case 2: e = integral.getX(Npc - 1 - j); n = 1;
                    break;
                case 3: e = -1.; n = integral.getX(Npc - 1 - j);
                    break;
                default:
                    cerr<<"Cos poszlo nie tak w obliczaniu warunkow brzegowych"<<endl;
                    exit(-1);
            }

           surface[k].operator()(j, 0) = N1(e, n);
           surface[k].operator()(j, 1) = N2(e, n);
           surface[k].operator()(j, 2) = N3(e, n);
           surface[k].operator()(j, 3) = N4(e, n);
        }
    }
}

void UniElement::calculateH_C(Element& e, const GlobalData& gd) {

    //zmienne lokalne - wektor macierzy H, dx i dy
    matrix dX(1, dN);
    matrix dY(1, dN);
    matrix resultX, resultY, H_pc, H;   //H
    matrix C_pc, C;                     //C
    //vector <matrix> H_pc(node , matrix(dN, dN));
    //jacobian jako wektor bo nie wiem czy beda potrzebne osobno czy nie
    vector <Jacobian> jacobians(node, Jacobian(2, 2));

    //--------------------------OBLICZANIE JAKOBIANU ---------------------------
    for (int i = 0; i < node; ++i) {            //node = tyle ile jakobianow - wiersze
        for (int j = 0; j < dN; ++j) {           //4 bo jakobian 2x2 - kolumny
            jacobians[i].J(0, 0) += (ksi_dN(i, j) * e.nodes[j].x);
            jacobians[i].J(0, 1) += (ksi_dN(i, j) * e.nodes[j].y);
            jacobians[i].J(1, 0) += (eta_dN(i, j) * e.nodes[j].x);
            jacobians[i].J(1, 1) += (eta_dN(i, j) * e.nodes[j].y);
        }
        jacobians[i].detJ = det(jacobians[i].J);
        jacobians[i].Jinv = inv(jacobians[i].J);
//        cout<<endl<<"Jakobian["<<i<<"]"<<jacobians[i].J<<endl
//        <<"Wyznacznik["<<i<<"] = "<<jacobians[i].detJ<<endl;
    }

    //--------------------------OBLICZANIE dN/dx i dN/dy oraz H oraz C---------------------------
    int N = -1;     //numer wiersza
    int M = -1;     //numer kolumny
    for (int i = 0; i < node; ++i) {

        if(N == (Npc - 1)) N = 0;
        else N++;
        if(i%Npc == 0) M++;

        for (int j = 0; j < dN; ++j) {
            dX(0, j) = jacobians[i].Jinv(0, 0) * ksi_dN(i, j) + jacobians[i].Jinv(0, 1) * eta_dN(i, j);
            dY(0, j) = jacobians[i].Jinv(1, 0) * ksi_dN(i, j) + jacobians[i].Jinv(1, 1) * eta_dN(i, j);
        }

        resultX = trans(get_row(dX,0)) * get_row(dX,0);
        resultY = trans(get_row(dY,0)) * get_row(dY,0);

        H_pc = gd.Conductivity * jacobians[i].detJ * (resultX + resultY);     //[H] = k(t) * ([dx] + [dy]) * dV
        H = H + (integral.getA(N) * integral.getA(M)) * H_pc;          // ostateczna macierz H

        //--------------- C ----------------
        //wszystko co powyzej dotyczy macierzy H!!
        C_pc = gd.Density * gd.SpecificHeat * jacobians[i].detJ * (trans(get_row(ksiEta_N,i)) * get_row(ksiEta_N, i));
        C = C + (integral.getA(N) * integral.getA(M)) * C_pc;
    }
//    cout<<"Macierz H["<<e.element_id<<"]"<<endl<<H<<endl;
    e.H = H;
    e.C = C;
}

void UniElement::calculateHbc_P(Element& e, const GlobalData& gd){

    matrix tmp, H_pc, Hbc;
    matrix P(dN, 1), P_pc;
    int first, second;

    for (int i = 0; i < dN; ++i) {
        //ustawianie odpowiednich punktów
        first = i;
        if(i == (dN-1)){
            second = 0;
        }else{
            second = i+1;
        }

        //Tylko dla brzegowych krawedzi
        if((e.nodes[first].BC == 1) && (e.nodes[second].BC == 1)){
            double detJ = sqrt(pow((e.nodes[first].x - e.nodes[second].x), 2) +
                               pow((e.nodes[first].y - e.nodes[second].y), 2)) / 2.;        // L / 2
            //cout<<"detJ: "<<detJ<<endl;

            for (int j = 0; j < Npc; ++j) {

                //Hbc
                tmp = trans(get_row(surface[i],j)) * get_row(surface[i],j);
                H_pc = H_pc + integral.getA(j) * tmp;
                //cout<<endl<<"i: "<<i<<endl<<"HPC LOKAL"<<endl<<H_pc;

                //P
                P_pc = P_pc + integral.getA(j) * trans(get_row(surface[i],j));
            }

            double mult = gd.Alfa * detJ;
            Hbc =  mult * H_pc + Hbc;
            P = gd.Tot *  mult * P_pc + P;
            //cout<<endl<<"HBC LOKAL"<<endl<<Hbc;
            //cout<<endl<<"P LOKAL"<<endl<<P;
            H_pc = H_pc - H_pc;                     //zerowanie
            P_pc = P_pc - P_pc;
        }
    }
    e.P = P;
    e.Hbc = Hbc;
}

void UniElement::printUniElements() const{
    cout<<endl<<"Macierz Ksi_dN (dla H):"<<endl;
    cout << ksi_dN << endl;
    cout<<"\nMacierz Eta_dN (dla H):"<<endl;
    cout << eta_dN << endl;
    cout<<"\nMacierz KsiEta_N (dla C):"<<endl;
    cout << ksiEta_N << endl;
    for (auto & s : surface) {
            cout<<endl<<"Macierz nodesN:"<<s<<endl;
    }
}

