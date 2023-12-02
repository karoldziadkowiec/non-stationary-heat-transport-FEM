#ifndef UNIVERSALELEMENT_H
#define UNIVERSALELEMENT_H
#include <iostream>
#include <cmath>

#include "GaussQuadrature.h"

using namespace std;

struct Surface
{
    double** ksiEtaMatrix;
    double** Ni;
    double*** Hbc_AtPci;
};

struct UniversalElement
{
    int N;
    double** dN_dKsi;
    double** dN_dEta;
    
    Surface surface[4];

    double** Ni_MatrixC;

    UniversalElement(int N);
    ~UniversalElement();

    void calculateShapeFunctionDerivatives();
    void printShapeFunctionDerivatives();

    void calculateKsiEtaMatrix_Values();
    void printKsiEtaMatrix_Values();
    void calculateMatrixOfN_Values(int surf);
    void printMatrixOfN_Values(int surf);


};

double dN1_dKsi(double eta);
double dN2_dKsi(double eta);
double dN3_dKsi(double eta);
double dN4_dKsi(double eta);
double dN1_dEta(double ksi);
double dN2_dEta(double ksi);
double dN3_dEta(double ksi);
double dN4_dEta(double ksi);

double N_Function(double ksi, double eta, int i);

#endif
