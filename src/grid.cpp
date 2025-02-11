//
// Created by Jakub Nowak on 18.10.2023.
//
#include "../include/grid.h"

//-------------------------------------GlobaData----------------------------------
void GlobalData::printGlobalData() const{
    cout<<"-----------------------------------GLOBAL DATA----------------------------------"<<endl<<
        "SimulationTime: "<<this->SimulationTime<<endl<<
        "SimulationStepTime: "<<this->SimulationStepTime<<endl<<
        "Conductivity: "<<this->Conductivity<<endl<<
        "Alfa: "<<this->Alfa<<endl<<
        "Tot: "<<this->Tot<<endl<<
        "InitialTemp: "<<this->InitialTemp<<endl<<
        "Density: "<<this->Density<<endl<<
        "SpecificHeat: "<<this->SpecificHeat<<endl<<
        "Nodes number: "<<this->Nodes_number<<endl<<
        "Elements number: "<<this->Elements_number
        <<endl;
}
GlobalData::GlobalData(){
    SimulationTime = SimulationStepTime = Conductivity =
    Alfa = Tot = InitialTemp = Density = SpecificHeat =
    Nodes_number = Elements_number = 0;
    }

//-------------------------------------Node---------------------------------------
Node::Node(){
    x = 0.;
    y = 0.;
    node_id = 0;
}

Node::Node(int ID, double wsp_x, double wsp_y){
    x = wsp_x;
    y = wsp_y;
    node_id = ID;
}

ostream& operator<<(ostream& output, const Node& obiekt) {
    // Implementacja wyjścia do strumienia
    output <<"("<< obiekt.node_id << ") ["<<setw(12)<<setprecision(9) <<
           obiekt.x << ",\t"<< setw(8) <<obiekt.y <<"]" <<"\tBC: "<<obiekt.BC<<endl;
    return output;
}

//-------------------------------------Element------------------------------------
Element::Element(): H(), Hbc(), C(), P(node_number, 1){
    element_id = 0;
    nodes[0] = Node();
    nodes[1] = Node();
    nodes[2] = Node();
    nodes[3] = Node();
}

Element::Element(int ID, Node ur, Node ul, Node dl, Node dr): H(), Hbc(), C(), P(node_number, 1){
    element_id = ID;
    nodes[0] = ur;
    nodes[1] = ul;
    nodes[2] = dl;
    nodes[3] = dr;
}

// dla samego poczaku -> czas = 0
void Element::calculateGlobalH(const GlobalData& gd) {
    H = H + Hbc;
}

//matrix t0 = matrix(node_number, 1, gd.InitialTemp);      // wektor pionowy t0
//matrix tmp = C/gd.SimulationStepTime;                   // macierz tymczasowa
//H = H + Hbc + tmp;                                     // ([H] * [C]/dT) * {t1} - układ rownan      //T=tał(czas)
//P = tmp * t0 + P;                                     // ([C]/dT)*{t0} + {P} - wyrazy wolne

ostream& operator<<(ostream& output, const Element& obiekt) {
    // Implementacja wyjścia do strumienia
    output <<endl<< "Element: " << obiekt.element_id << endl <<
           obiekt.nodes[0] << obiekt.nodes[1] <<
           obiekt.nodes[2] << obiekt.nodes[3] <<endl
           <<"Macierz H:"<<obiekt.H<<endl<<
           //endl<<"Macierz Hbc:"<<obiekt.Hbc<<endl<<
           endl<<"Wektor P:"<<obiekt.P<<endl<<
            endl<<"Macierz C:"<<obiekt.C<<endl;

    return output;
}

//-------------------------------------Grid---------------------------------------
void Grid::printGrid(){
    cout<<endl<<endl<<"--------------------------------------GRID-------------------------------------"<<endl;
    for (const Element& element : this->elements) {
        std::cout << element << " ";
        cout<<"----------------------------------------------------------"<<endl;
    }
}
//---------------------------------wczytywanie z pliku-------------------------------
void readFromFile(ifstream& file, GlobalData& gd, Grid& grid){

    //deklaracje struktur tymczasowych
    Node temp_node;
    Element temp_element;

    if (file.is_open()) {

        //potrzebne zmienne
        string line;
        string uselessLine;
        string key;
        int value;
        char coma;
        int ID,A, B, C, D;
        short flagBC;
        double X, Y;


        //wczytywanie danych globalnych------------------------------------
        int counter=1;
        for (int i=0;i<10;i++) {

            getline(file, line);
            istringstream stream(line);
            stream >> key >> value;

            switch(counter){
                case 1: gd.SimulationTime = value; break;
                case 2: gd.SimulationStepTime = value; break;
                case 3: gd.Conductivity = value; break;
                case 4: gd.Alfa = value; break;
                case 5: gd.Tot = value; break;
                case 6: gd.InitialTemp = value; break;
                case 7: gd.Density = value; break;
                case 8: gd.SpecificHeat = value; break;
                case 9: gd.Nodes_number = value; break;
                case 10: gd.Elements_number = value; break;
                default: cout <<" Koniec danych globalnych! "<<endl;
            }
            counter++;
        }

        //wczytywanie nodes-------------------------------------------
        getline(file, uselessLine); //wczytanie linii *Nodes

        for (int i = 0; i < gd.Nodes_number; ++i){

            getline(file, line);
            istringstream stream(line);
            stream >> ID >> coma >> X >> coma >> Y;

            temp_node.node_id = ID;
            temp_node.x = X;
            temp_node.y = Y;

            grid.nodes.push_back(temp_node);
            //temp_node.print_node();
        }

        //wczytywanie elements-----------------------------------------
        getline(file, uselessLine); //wczytanie linii *Elements

        for (int i = 0; i < gd.Elements_number; ++i){
            getline(file, line);
            istringstream stream(line);
            stream >> ID >> coma >> A >> coma >> B >> coma >> C >> coma >> D;

            temp_element.element_id = ID;
            temp_element.nodes[0].node_id = A;
            temp_element.nodes[1].node_id = B;
            temp_element.nodes[2].node_id = C;
            temp_element.nodes[3].node_id = D;

//            temp_element.nodes[0] = grid.nodes[A-1];
//            temp_element.nodes[1] = grid.nodes[B-1];
//            temp_element.nodes[2] = grid.nodes[C-1];
//            temp_element.nodes[3] = grid.nodes[D-1];

            grid.elements.push_back(temp_element);
        }

        //wczytywanie flag BC-----------------------------------------
        getline(file, uselessLine); //wczytanie linii *BC


        while (getline(file, line, ',')) {
            istringstream stream(line);
            stream >> flagBC;
            grid.nodes[flagBC-1].BC = 1;                  //rozwiazanie tymczasowe
        }

        //przypisywanie nodow do elements
        for (int i = 0; i < gd.Elements_number; ++i){
            for(auto & node : grid.elements[i].nodes){
                node = grid.nodes[node.node_id-1];
            }
        }

    }else{
        cerr<<"Plik jest zamkniety"<<endl;
    }
}

//zapis do pliku przeznaczony dla programu ParaView
void writeToFile(ofstream& file, const GlobalData& gd, const Grid& grid, const matrix& X){
    int elType = 9;                     //typ elementu dla paraview
    file << "# vtk DataFile Version 2.0\n";
    file << "Unstructured Grid Example\n";
    file << "ASCII\n";
    file << "DATASET UNSTRUCTURED_GRID\n\n";

    file << "POINTS "<< gd.Nodes_number <<" float\n";
    for (int i = 0; i < gd.Nodes_number; ++i) {
        file<<grid.nodes[i].x<<" "<<grid.nodes[i].y<<" "<<0<<"\n";
    }

    file << "\nCELLS "<<gd.Elements_number<<" "<<gd.Elements_number*5<<"\n";     //chyba * 5
    for (int i = 0; i < gd.Elements_number; ++i) {
        file<<4;
        for (int j = 0; j < 4; ++j) {
            file<<" "<<grid.elements[i].nodes[j].node_id-1;
        }
        file<<"\n";
    }

    file<<"\nCELL_TYPES "<<gd.Elements_number<<"\n";
    for (int i = 0; i < gd.Elements_number; ++i) {
        file <<elType<<"\n";
    }

    file << "\nPOINT_DATA "<<gd.Nodes_number<<"\n";
    file << "SCALARS Temp float 1\n";
    file << "LOOKUP_TABLE default\n";
    for (int i = 0; i < gd.Nodes_number; ++i) {
       file << X(i,0) <<"\n";
    }
}
