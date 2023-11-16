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

    double*** xH_AtPci;
    double*** yH_AtPci;

    double*** Hpci;

    double** H;

    Jakobian(int N);
    ~Jakobian();

    void calculateDerivativesAtPci(const UniversalElement& universalElement, const Grid& grid, int elementNumber, int pc);
    void printJakobianMatrix();
    double calculateDetJ();
    void printDetJ();
    double calculate1_DetJ();
    void calculateJakobianMatrix();
    void calculateShapeFunctionDerivativesForPci(const UniversalElement& universalElement, int pc);
    void printShapeFunctionDerivativesForPci(int pc);
    void calculateMatrixHForXandYForPci();
    void printMatrixHForXandYForPci();
    void calculateMatrixHpci(int conductivity);
    void printMatrixHpci();
    void calculateMatrixH();
    void printMatrixH();
};

#endif