#ifndef JAKOBIAN_H
#define JAKOBIAN_H
#include <iostream>
#include <cmath>

#include "GaussQuadrature.h"
#include "UniversalElement.h"

using namespace std;

struct Jakobian
{
    // Elementy macierzy Jakobiego
    double dx_dKsi;
    double dy_dKsi;
    double dx_dEta;
    double dy_dEta;

    double** dN_dx;
    double** dN_dy;
    int N;

    double kt; // Wsp�czynnik przewodno�ci cieplnej
    double dV; // Obj�to�� elementu
    double dS; // Powierzchnia elementu

    double*** xH_AtPci;
    double*** yH_AtPci;

    double*** Hpci;

    double*** Hbc_AtPci;
    double*** Hbci;

    double*** Cpci;

    Jakobian(int N);
    ~Jakobian();

    //H
    void calculateDerivativesAtPci(const UniversalElement& universalElement, const Grid& grid, int elementNumber, int pc);
    void printJakobianMatrix();
    double calculateDetJ();
    void printDetJ();
    double calculate1_DetJ();
    void calculateJakobianMatrix();
    void calculateShapeFunctionDerivativesForPci(const UniversalElement& universalElement, int pc);
    void printShapeFunctionDerivativesForPci(int pc);
    void calculateMatrixHForXandYForPci(int pc);
    void calculateMatrixHpci(int pc, int conductivity);
    void printMatrixHpci();
    void calculateMatrixH(const Grid& grid, int elementNumber);
    void printMatrixH(const Grid& grid, int elementNumber);

    //Hbc
    double calculateHbcDetJ(const Grid& grid, int Nx, int Nk);
    void printHbcDetJ(const Grid& grid, int Nx, int Nk);
    void calculateMatrixHbciForPci(const UniversalElement& universalElement, int surface);
    void calculateMatrixHbci(int surface, int alfa, const Grid& grid, int Nx, int Nk);
    void printMatrixHbci();
    void calculateMatrixHbc(const Grid& grid, int elementNumber);
    void printMatrixHbc(const Grid& grid, int elementNumber);
    void zeroMatrixHbci();
    void sumMatrixH_Hbc(const Grid& grid, int elementNumber);

    //Vector P
    void zeroVectorP(const Grid& grid, int elementNumber);
    void calculateVectorP_ForPci(const Grid& grid, const UniversalElement& universalElement, int surface, int elementNumber, int tot, int alfa, int Nx, int Nk);
    void printVectorP(const Grid& grid, int elementNumber);

    //C
    void calculateMatrixCpci(const UniversalElement& universalElement, int pc, int specificHeat, int density);
    void printMatrixCpci();
    void calculateMatrixC(const Grid& grid, int elementNumber);
    void printMatrixC(const Grid& grid, int elementNumber);
};

#endif