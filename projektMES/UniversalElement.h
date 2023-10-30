#ifndef UNIVERSALELEMENT_H
#define UNIVERSALELEMENT_H
#include <iostream>
#include <cmath>

#include "GaussQuadrature.h"

using namespace std;

struct UniversalElement
{
    int N;
    double** dN_dKsi;
    double** dN_dEta;

    UniversalElement(int N);
    ~UniversalElement();

    void computeShapeFunctionDerivatives();
    void printShapeFunctionDerivatives();
};

double dN1_dKsi(double eta);
double dN2_dKsi(double eta);
double dN3_dKsi(double eta);
double dN4_dKsi(double eta);
double dN1_dEta(double ksi);
double dN2_dEta(double ksi);
double dN3_dEta(double ksi);
double dN4_dEta(double ksi);

#endif