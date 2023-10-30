#include <iostream>
#include <fstream>

#include "GlobalData.h"
#include "Grid.h"
#include "GaussQuadrature.h"
#include "UniversalElement.h"

using namespace std;

//LAB 2
double function_1D(double x);
double function_2D(double x, double y);
void lab2();
//LAB3
void lab3();

int main()
{
    Grid grid;
    GlobalData globalData;
    string fileName = "Test1_4_4.txt";

    readDataFromFile(fileName, globalData, grid);
    printGridData(globalData, grid);

    lab2();
    lab3();

    delete[] grid.nodes;
    delete[] grid.elements;

    return 0;
}

//LAB 2
double function_1D(double x) {
    return 5 * pow(x, 2) + 3 * x + 6;
}

double function_2D(double x, double y) {
    return 5 * pow(x, 2) * pow(y, 2) + 3 * x * y + 6;
}

void lab2() {
    int N; //Nodes number

    cout << "\nIntegration interval: [-1,1]" << endl;
    cout << "f(x) = 5 * x^2 + 3 * x + 6" << endl;
    //Gaussian integration in 1D space using 2 points integration scheme
    N = 2;
    double integral = gauss(N, function_1D);
    cout << "Result of 2-point Gauss Quadrature in 1D: " << integral << endl;

    //Gaussian integration in 1D space using 3 points integration scheme
    N = 3;
    double integral2 = gauss(N, function_1D);
    cout << "Result of 3-point Gauss Quadrature in 1D: " << integral2 << endl << endl;

    cout << "f(x) = 5 * x^2 * y^2 + 3 * x * y + 6" << endl;
    //Gaussian integration in 2D space using 2 points integration scheme
    N = 2;
    double integral3 = gauss(N, function_2D);
    cout << "Result of 2-point Gauss Quadrature in 2D: " << integral3 << endl;

    //Gaussian integration in 2D space using 3 points integration scheme
    N = 3;
    double integral4 = gauss(N, function_2D);
    cout << "Result of 3-point Gauss Quadrature in 2D: " << integral4 << endl;
}

//LAB3
void lab3() {
    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.computeShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();
}