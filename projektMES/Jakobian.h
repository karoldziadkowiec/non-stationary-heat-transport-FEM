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

    double** xH_AtPc1;
    double** yH_AtPc1;
    double** xH_AtPc2;
    double** yH_AtPc2;
    double** xH_AtPc3;
    double** yH_AtPc3;
    double** xH_AtPc4;
    double** yH_AtPc4;

    double** Hpc1;
    double** Hpc2;
    double** Hpc3;
    double** Hpc4;

    double** H;

    Jakobian(int N);
    ~Jakobian();

    void calculateDerivativesAtPci(const UniversalElement& universalElement, const Grid& grid, int i);
    void printJakobianMatrix();
    double calculateDetJ();
    double calculate1_DetJ();
    void calculateJakobianMatrix();
    void calculateShapeFunctionDerivativesForPci(const UniversalElement& universalElement);
    void printShapeFunctionDerivativesForPci();
    void calculateMatrixHForXandYForPci();
    void printMatrixHForXandYForPci();
    void calculateMatrixHpci(int kt);
    void printMatrixHpci();
    void calculateMatrixH();
    void printMatrixH();
};

#endif