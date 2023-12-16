#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>
#include <sys/stat.h> // Dodane dla funkcji mkdir w systemach Unix
#include <direct.h>   // Dodane dla funkcji _mkdir w systemie Windows

#include "ParaView.h"

using namespace std;

bool createDirectory(const string& path) // Tworzenie katalogu dla konkretnej siatki
{
#ifdef _WIN32 // sprawdza, czy kod jest kompilowany pod systemem Windows
    return _mkdir(path.c_str()) == 0; // Funkcja _mkdir zwraca 0, jeœli operacja utworzenia folderu siê powiedzie
#else // wersja Unix
    return mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0; //Opcje S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH ustawiaj¹ prawa dostêpu do folderu, umo¿liwiaj¹c pe³ny dostêp dla w³aœciciela i odczyt/uruchamianie dla innych. 
#endif
}

// Tworzenie pliku ParaView dla konkretnej siatki
void createParaViewFile(string katalogName, GlobalData& globalData, Grid& grid, SoE& soe, int time)
{
    createDirectory(katalogName);
    string filename = katalogName + "/foo" + to_string(time) + ".vtk";
    ofstream file(filename);

    if (!file.is_open())
    {
        cerr << "The file cannot be opened." << endl;
        return;
    }

    file << "# vtk DataFile Version 2.0" << endl;
    file << "Unstructured Grid Example" << endl;
    file << "ASCII" << endl;
    file << "DATASET UNSTRUCTURED_GRID" << endl<<endl;

    file << "POINTS "<< globalData.nodesNumber <<" float" << endl;
    for (int i = 0; i < globalData.nodesNumber; i++)
    {
        file << grid.nodes[i].x << " " << grid.nodes[i].y << " 0" <<endl;
    }

    file << "\nCELLS " << globalData.elementsNumber << " " << globalData.elementsNumber * 5 << endl;
    for (int i = 0; i < globalData.elementsNumber; i++)
    {
        file << "4 " << grid.elements[i].id[0] - 1 << " " << grid.elements[i].id[1] - 1 << " " << grid.elements[i].id[2] - 1 << " " << grid.elements[i].id[3] - 1 << endl;
    }

    file << "\nCELL_TYPES " << globalData.elementsNumber << endl;
    for (int i = 0; i < globalData.elementsNumber; i++)
    {
        file << "9" << endl;
    }

    file << "\nPOINT_DATA " << globalData.nodesNumber << endl;
    file << "SCALARS Temp float 1" << endl;
    file << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < globalData.nodesNumber; i++)
    {
        file << soe.t[i] << endl;
    }

    file.close();
}