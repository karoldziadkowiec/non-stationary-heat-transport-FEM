#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H

using namespace std;

struct GaussQuadrature
{
    int* k; // gaussian point index
    double* xk; // coordinate of the integration point
    double* wk; // weight of the integration point
};

GaussQuadrature returnRowOfGaussTable(int N);
double gauss(int N, double (*function)(double));
double gauss(int N, double (*function)(double, double));

#endif
