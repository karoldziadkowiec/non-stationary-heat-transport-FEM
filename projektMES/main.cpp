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
//LAB5
void lab5();
void lab5_Test1_4_4();
void lab5_Test2_4_4_MixGrid();

int main()
{
    //lab1();
    //lab2();
    //lab3();
    //lab4();
    //lab4_Test1_4_4();
    //lab4_Test2_4_4_MixGrid();
    lab5();
    //lab5_Test1_4_4();
    //lab5_Test2_4_4_MixGrid();

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
    int kt = 30; // conductivity

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
    universalElement.printShapeFunctionDerivatives();

    Jakobian jakobian(N);
    for (int elNumber = 0; elNumber < globalData.elementsNumber; elNumber++) {
        cout << "\n\n\t\tELEMENT " << elNumber + 1 << endl;
        for (int pc = 0; pc < N * N; pc++) {
            cout << "\n\tPunkt calkowania " << pc + 1 << endl;
            jakobian.calculateDerivativesAtPci(universalElement, testGrid, elNumber, pc);
            jakobian.printJakobianMatrix();
            jakobian.printDetJ();
            jakobian.calculateJakobianMatrix();
            jakobian.calculateShapeFunctionDerivativesForPci(universalElement, pc);
            jakobian.printShapeFunctionDerivativesForPci(pc);
            jakobian.calculateMatrixHForXandYForPci(pc);
            jakobian.calculateMatrixHpci(pc, kt);
        }
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH(testGrid, elNumber);
        jakobian.printMatrixH(testGrid, elNumber);
    }
}

void lab4_Test1_4_4() {
    cout << "\n Test1_4_4.txt" << endl;

    Grid grid;
    GlobalData globalData;
    string fileName = "Test1_4_4.txt";

    readDataFromFile(fileName, globalData, grid);
    printGridData(globalData, grid);
    int kt = globalData.conductivity; // conductivity

    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.calculateShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();

    Jakobian jakobian(N);
    for (int elNumber = 0; elNumber < globalData.elementsNumber; elNumber++) {
        cout << "\n\n\t\tELEMENT " << elNumber + 1 << endl;
        for (int pc = 0; pc < N * N; pc++) {
            cout << "\n\tPunkt calkowania " << pc + 1 << endl;
            jakobian.calculateDerivativesAtPci(universalElement, grid, elNumber, pc);
            jakobian.printJakobianMatrix();
            jakobian.printDetJ();
            jakobian.calculateJakobianMatrix();
            jakobian.calculateShapeFunctionDerivativesForPci(universalElement, pc);
            jakobian.printShapeFunctionDerivativesForPci(pc);
            jakobian.calculateMatrixHForXandYForPci(pc);
            jakobian.calculateMatrixHpci(pc, kt);
        }
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH(grid, elNumber);
        jakobian.printMatrixH(grid, elNumber);
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
    int kt = globalData.conductivity; // conductivity

    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.calculateShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();

    Jakobian jakobian(N);
    for (int elNumber = 0; elNumber < globalData.elementsNumber; elNumber++) {
        cout << "\n\n\t\tELEMENT " << elNumber + 1 << endl;
        for (int pc = 0; pc < N * N; pc++) {
            cout << "\n\tPunkt calkowania " << pc + 1 << endl;
            jakobian.calculateDerivativesAtPci(universalElement, grid, elNumber, pc);
            jakobian.printJakobianMatrix();
            jakobian.printDetJ();
            jakobian.calculateJakobianMatrix();
            jakobian.calculateShapeFunctionDerivativesForPci(universalElement, pc);
            jakobian.printShapeFunctionDerivativesForPci(pc);
            jakobian.calculateMatrixHForXandYForPci(pc);
            jakobian.calculateMatrixHpci(pc, kt);
        }
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH(grid, elNumber);
        jakobian.printMatrixH(grid, elNumber);
    }

    delete[] grid.nodes;
    delete[] grid.elements;
}

void lab5() {
    cout << "\nTEST" << endl;
    Grid testGrid;
    GlobalData globalData;
    globalData.elementsNumber = 1;
    int kt = 30; // conductivity
    int alfa = 25; // heat transfer coefficient

    testGrid.elements = new Element[1];
    testGrid.elements[0].id[0] = 1;
    testGrid.elements[0].id[1] = 2;
    testGrid.elements[0].id[2] = 4;
    testGrid.elements[0].id[3] = 3;

    testGrid.nodes = new Node[4];
    testGrid.nodes[3].x = 0;
    testGrid.nodes[3].y = 0;

    testGrid.nodes[2].x = 0.025;
    testGrid.nodes[2].y = 0;

    testGrid.nodes[0].x = 0.025;
    testGrid.nodes[0].y = 0.025;

    testGrid.nodes[1].x = 0;
    testGrid.nodes[1].y = 0.025;

    testGrid.nodes[0].BC = 1;
    testGrid.nodes[1].BC = 1;
    testGrid.nodes[2].BC = 1;
    testGrid.nodes[3].BC = 1;

    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.calculateShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();

    universalElement.calculateKsiEtaMatrix_Values();
    universalElement.printKsiEtaMatrix_Values();

    Jakobian jakobian(N);
    for (int elNumber = 0; elNumber < globalData.elementsNumber; elNumber++) {
        cout << "\n\n\t\tELEMENT " << elNumber + 1 << endl;
        for (int pc = 0; pc < N * N; pc++) {
            cout << "\n\tPunkt calkowania " << pc + 1 << endl;
            jakobian.calculateDerivativesAtPci(universalElement, testGrid, elNumber, pc);
            jakobian.printJakobianMatrix();
            jakobian.printDetJ();
            jakobian.calculateJakobianMatrix();
            jakobian.calculateShapeFunctionDerivativesForPci(universalElement, pc);
            jakobian.printShapeFunctionDerivativesForPci(pc);
            jakobian.calculateMatrixHForXandYForPci(pc);
            jakobian.calculateMatrixHpci(pc, kt);
        }
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH(testGrid, elNumber);
        jakobian.printMatrixH(testGrid, elNumber);

        for (int surface = 0; surface < 4; surface++) {
            if (surface == 0) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[2] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[3] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                }
            }
            else if (surface == 1) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[3] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[0] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                }
            }
            else if (surface == 2) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[0] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[1] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                }
            }
            else if (surface == 3) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[1] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[2] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                }
            }
        }
        jakobian.printMatrixHbci();
        jakobian.calculateMatrixHbc(testGrid, elNumber);
        jakobian.printMatrixHbc(testGrid, elNumber);
        jakobian.zeroMatrixHbci();
    }
}

void lab5_Test1_4_4() {
    cout << "\n Test1_4_4.txt" << endl;

    Grid testGrid;
    GlobalData globalData;
    string fileName = "Test1_4_4.txt";

    readDataFromFile(fileName, globalData, testGrid);
    printGridData(globalData, testGrid);
    int kt = globalData.conductivity; // conductivity
    int alfa = globalData.alfa; // heat transfer coefficient

    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.calculateShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();

    universalElement.calculateKsiEtaMatrix_Values();
    universalElement.printKsiEtaMatrix_Values();

    Jakobian jakobian(N);
    for (int elNumber = 0; elNumber < globalData.elementsNumber; elNumber++) {
        cout << "\n\n\t\tELEMENT " << elNumber + 1 << endl;
        for (int pc = 0; pc < N * N; pc++) {
            cout << "\n\tPunkt calkowania " << pc + 1 << endl;
            jakobian.calculateDerivativesAtPci(universalElement, testGrid, elNumber, pc);
            jakobian.printJakobianMatrix();
            jakobian.printDetJ();
            jakobian.calculateJakobianMatrix();
            jakobian.calculateShapeFunctionDerivativesForPci(universalElement, pc);
            jakobian.printShapeFunctionDerivativesForPci(pc);
            jakobian.calculateMatrixHForXandYForPci(pc);
            jakobian.calculateMatrixHpci(pc, kt);
        }
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH(testGrid, elNumber);
        jakobian.printMatrixH(testGrid, elNumber);

        for (int surface = 0; surface < 4; surface++) {
            if (surface == 0) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[2] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[3] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                }
            }
            else if (surface == 1) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[3] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[0] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                }
            }
            else if (surface == 2) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[0] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[1] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                }
            }
            else if (surface == 3) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[1] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[2] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                }
            }
        }
        jakobian.printMatrixHbci();
        jakobian.calculateMatrixHbc(testGrid, elNumber);
        jakobian.printMatrixHbc(testGrid, elNumber);
        jakobian.zeroMatrixHbci();
    }
}

void lab5_Test2_4_4_MixGrid()
{
    cout << "\n Test2_4_4_MixGrid.txt" << endl;

    Grid testGrid;
    GlobalData globalData;
    string fileName = "Test2_4_4_MixGrid.txt";

    readDataFromFile(fileName, globalData, testGrid);
    printGridData(globalData, testGrid);
    int kt = globalData.conductivity; // conductivity
    int alfa = globalData.alfa; // heat transfer coefficient

    int N = 2; //Nodes number
    UniversalElement universalElement(N);
    universalElement.calculateShapeFunctionDerivatives();
    universalElement.printShapeFunctionDerivatives();

    universalElement.calculateKsiEtaMatrix_Values();
    universalElement.printKsiEtaMatrix_Values();

    Jakobian jakobian(N);
    for (int elNumber = 0; elNumber < globalData.elementsNumber; elNumber++) {
        cout << "\n\n\t\tELEMENT " << elNumber + 1 << endl;
        for (int pc = 0; pc < N * N; pc++) {
            cout << "\n\tPunkt calkowania " << pc + 1 << endl;
            jakobian.calculateDerivativesAtPci(universalElement, testGrid, elNumber, pc);
            jakobian.printJakobianMatrix();
            jakobian.printDetJ();
            jakobian.calculateJakobianMatrix();
            jakobian.calculateShapeFunctionDerivativesForPci(universalElement, pc);
            jakobian.printShapeFunctionDerivativesForPci(pc);
            jakobian.calculateMatrixHForXandYForPci(pc);
            jakobian.calculateMatrixHpci(pc, kt);
        }
        jakobian.printMatrixHpci();
        jakobian.calculateMatrixH(testGrid, elNumber);
        jakobian.printMatrixH(testGrid, elNumber);

        for (int surface = 0; surface < 4; surface++) {
            if (surface == 0) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[2] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[3] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[2] - 1, testGrid.elements[elNumber].id[3] - 1);
                }
            }
            else if (surface == 1) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[3] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[0] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[3] - 1, testGrid.elements[elNumber].id[0] - 1);
                }
            }
            else if (surface == 2) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[0] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[1] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[0] - 1, testGrid.elements[elNumber].id[1] - 1);
                }
            }
            else if (surface == 3) {
                if ((testGrid.nodes[testGrid.elements[elNumber].id[1] - 1].BC != 0) && (testGrid.nodes[testGrid.elements[elNumber].id[2] - 1].BC != 0)) {
                    universalElement.calculateMatrixOfN_Values(surface);
                    universalElement.printMatrixOfN_Values(surface);
                    jakobian.calculateHbcDetJ(testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                    jakobian.printHbcDetJ(testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                    jakobian.calculateMatrixHbciForPci(universalElement, surface);
                    jakobian.calculateMatrixHbci(surface, alfa, testGrid, testGrid.elements[elNumber].id[1] - 1, testGrid.elements[elNumber].id[2] - 1);
                }
            }
        }
        jakobian.printMatrixHbci();
        jakobian.calculateMatrixHbc(testGrid, elNumber);
        jakobian.printMatrixHbc(testGrid, elNumber);
        jakobian.zeroMatrixHbci();
    }
}