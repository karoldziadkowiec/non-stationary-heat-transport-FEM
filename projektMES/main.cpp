#include <iostream>
#include <fstream>

#include "GlobalData.h"
#include "Grid.h"
#include "GaussQuadrature.h"
#include "UniversalElement.h"
#include "Jakobian.h"

using namespace std;

//LAB1
void lab1();
//LAB 2
double function_1D(double x);
double function_2D(double x, double y);
void lab2();
//LAB3
void lab3();
//LAB4
void lab4();
void lab4_Test1_4_4();
void lab4_Test2_4_4_MixGrid();

int main()
{
    //lab1();
    //lab2();
    //lab3();
    lab4();
    lab4_Test1_4_4();
    //lab4_Test2_4_4_MixGrid();

    return 0;
}

//LAB 1
void lab1() {
    Grid grid;
    GlobalData globalData;
    string fileName = "Test1_4_4.txt";

    readDataFromFile(fileName, globalData, grid);
    printGridData(globalData, grid);

    delete[] grid.nodes;
    delete[] grid.elements;
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
    universalElement.calculateShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();
}

//LAB4
void lab4() {
    cout << "\nTEST" << endl;
    Grid testGrid;
    GlobalData globalData;
    globalData.elementsNumber = 1;
    int kt = 30;

    testGrid.elements = new Element[1];
    testGrid.elements[0].id[0] = 1;
    testGrid.elements[0].id[1] = 2;
    testGrid.elements[0].id[2] = 3;
    testGrid.elements[0].id[3] = 4;

    testGrid.nodes = new Node[4];
    testGrid.nodes[0].x = 0;
    testGrid.nodes[0].y = 0;

    testGrid.nodes[1].x = 0.025;
    testGrid.nodes[1].y = 0;

    testGrid.nodes[2].x = 0.025;
    testGrid.nodes[2].y = 0.025;

    testGrid.nodes[3].x = 0;
    testGrid.nodes[3].y = 0.025;

    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.calculateShapeFunctionDerivatives();

    Jakobian jakobian(N);
    for (int i = 0; i < globalData.elementsNumber; i++) {
        jakobian.calculateDerivativesAtPci(universalElement, testGrid, i);
        jakobian.printJakobianMatrix();
        jakobian.calculateJakobianMatrix();
        jakobian.printJakobianMatrix();
        jakobian.calculateShapeFunctionDerivativesForPci(universalElement);
        jakobian.printShapeFunctionDerivativesForPci();
        jakobian.calculateMatrixHForXandYForPci();
        jakobian.printMatrixHForXandYForPci();
        jakobian.calculateMatrixHpci(kt);
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH();
        jakobian.printMatrixH();
    }
}

void lab4_Test1_4_4() {
    cout << "\n Test1_4_4.txt" << endl;

    Grid grid;
    GlobalData globalData;
    string fileName = "Test1_4_4.txt";

    readDataFromFile(fileName, globalData, grid);
    printGridData(globalData, grid);
    int kt = globalData.conductivity;

    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.calculateShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();

    Jakobian jakobian(N);
    for (int i = 0; i < globalData.elementsNumber; i++) {
        cout << "\n\tELEMENT " << i + 1 << endl;
        jakobian.calculateDerivativesAtPci(universalElement, grid, i);
        jakobian.printJakobianMatrix();
        jakobian.calculateJakobianMatrix();
        jakobian.printJakobianMatrix();
        jakobian.calculateShapeFunctionDerivativesForPci(universalElement);
        jakobian.printShapeFunctionDerivativesForPci();
        jakobian.calculateMatrixHForXandYForPci();
        jakobian.printMatrixHForXandYForPci();
        jakobian.calculateMatrixHpci(kt);
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH();
        jakobian.printMatrixH();
    }

    delete[] grid.nodes;
    delete[] grid.elements;
}

void lab4_Test2_4_4_MixGrid() {
    cout << "\n Test2_4_4_MixGrid.txt" << endl;

    Grid grid;
    GlobalData globalData;
    string fileName = "Test2_4_4_MixGrid.txt";

    readDataFromFile(fileName, globalData, grid);
    printGridData(globalData, grid);
    int kt = globalData.conductivity;

    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.calculateShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();

    Jakobian jakobian(N);
    for (int i = 0; i < globalData.elementsNumber; i++) {
        cout << "\n\tELEMENT " << i + 1 << endl;
        jakobian.calculateDerivativesAtPci(universalElement, grid, i);
        jakobian.printJakobianMatrix();
        jakobian.calculateJakobianMatrix();
        jakobian.printJakobianMatrix();
        jakobian.calculateShapeFunctionDerivativesForPci(universalElement);
        jakobian.printShapeFunctionDerivativesForPci();
        jakobian.calculateMatrixHForXandYForPci();
        jakobian.printMatrixHForXandYForPci();
        jakobian.calculateMatrixHpci(kt);
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH();
        jakobian.printMatrixH();
    }

    delete[] grid.nodes;
    delete[] grid.elements;
}