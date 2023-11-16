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

    //Surface surface[4];

    UniversalElement(int N);
    ~UniversalElement();

    void calculateShapeFunctionDerivatives();
    void printShapeFunctionDerivatives();
};

struct Surface
{
    int n;
    double** N;

    Surface(int n);
    ~Surface();
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