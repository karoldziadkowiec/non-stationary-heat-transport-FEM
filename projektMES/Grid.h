#ifndef GRID_H
#define GRID_H

#include "GlobalData.h"
#include <string> 

struct Node
{
    double x;
    double y;
    int BC;
};

struct Element
{
    int id[4];
    double H[4][4];
    double Hbc[4][4];
    double P[4];
    double C[4][4];

    double HplusC_dT[4][4];
    double Ct0_dTplusP[4];
};

struct Grid
{
    Element* elements;
    Node* nodes;
};

void readDataFromFile(const std::string& filename, GlobalData& globalData, Grid& grid);
void printGridData(const GlobalData& globalData, const Grid& grid);

#endif
