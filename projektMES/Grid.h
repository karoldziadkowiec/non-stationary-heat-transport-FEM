#ifndef GRID_H
#define GRID_H

#include "GlobalData.h"

struct Node
{
    double x;
    double y;
};

struct Element
{
    int id[4];
};

struct Grid
{
    Element* elements;
    Node* nodes;
};

void readDataFromFile(const std::string& filename, GlobalData& globalData, Grid& grid);
void printGridData(const GlobalData& globalData, const Grid& grid);

#endif
