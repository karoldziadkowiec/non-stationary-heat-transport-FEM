#ifndef JAKOBIAN_H
#define JAKOBIAN_H
#include <iostream>
#include <cmath>

#include "GaussQuadrature.h"
#include "UniversalElement.h"

using namespace std;

struct Jakobian
{
    double dx_dKsi;
    double dy_dKsi;
    double dx_dEta;
    double dy_dEta;

    double** dN_dx;
    double** dN_dy;
    int N;

    double kt;
    double dV;
    double dS;

    double*** xH_AtPci;
    double*** yH_AtPci;

    double*** Hpci;

    double*** Hbc_AtPci;
    double*** Hbci;

    Jakobian(int N);
    ~Jakobian();

    //Hpc
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
};

#endif