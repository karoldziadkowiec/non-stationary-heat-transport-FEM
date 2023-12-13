#ifndef PARAVIEW_H
#define PARAVIEW_H
#include <iostream>
#include <cmath>
#include <filesystem>

#include "GlobalData.h"
#include "Grid.h"
#include "Agregation.h"

using namespace std;

bool createDirectory(const string& path);
void createParaViewFile(string katalogName, GlobalData& globalData, Grid& grid, SoE& soe, int time);

#endif