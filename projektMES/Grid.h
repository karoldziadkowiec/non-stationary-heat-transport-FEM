#ifndef GRID_H
#define GRID_H

#include "GlobalData.h"
#include <string> 

struct Node
{
    double x;
    double y;
    int bc;
};

struct Element
{
    int id[4];
    int H[4][4];
    int Hbc[4][4];
};

struct Grid
{
    Element* elements;
    Node* nodes;
};

void readDataFromFile(const std::string& filename, GlobalData& globalData, Grid& grid);
void printGridData(const GlobalData& globalData, const Grid& grid);

#endif
