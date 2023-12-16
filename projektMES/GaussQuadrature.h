#ifndef GAUSSQUADRATURE_H
#define GAUSSQUADRATURE_H

using namespace std;

struct GaussQuadrature
{
    int* k; // index wiersza
    double* xk; // punkt ca³kowania
    double* wk; // waga
};

GaussQuadrature returnRowOfGaussTable(int N);
double gauss(int N, double (*function)(double));
double gauss(int N, double (*function)(double, double));

#endif