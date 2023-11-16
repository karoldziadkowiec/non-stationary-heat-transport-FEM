#include <cmath>

#include "UniversalElement.h"

UniversalElement::UniversalElement(int N)
{
    this->N = N;

    dN_dKsi = new double* [N * N];
    dN_dEta = new double* [N * N];

    for (int i = 0; i < N * N; i++) {
        dN_dKsi[i] = new double[4] {};
        dN_dEta[i] = new double[4] {};
    }
}

UniversalElement::~UniversalElement()
{
    for (int i = 0; i < N * N; i++) {
        delete[] dN_dKsi[i];
        delete[] dN_dEta[i];
    }
    delete[] dN_dKsi;
    delete[] dN_dEta;
}

Surface::Surface(int n)
{
    this->n = n;

    N = new double* [n * n];

    for (int i = 0; i < n * n; i++) {
        N[i] = new double[4] {};
    }
}

Surface::~Surface()
{
    for (int i = 0; i < n * n; i++) {
        delete[] N[i];
    }
    delete[] N;
}

void UniversalElement::calculateShapeFunctionDerivatives()
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int pc = i * N + j;

            // dNi/dKsi
            dN_dKsi[pc][0] = dN1_dKsi(tableRow.xk[i]);
            dN_dKsi[pc][1] = dN2_dKsi(tableRow.xk[i]);
            dN_dKsi[pc][2] = dN3_dKsi(tableRow.xk[i]);
            dN_dKsi[pc][3] = dN4_dKsi(tableRow.xk[i]);

            // dNi/dEta
            dN_dEta[pc][0] = dN1_dEta(tableRow.xk[j]);
            dN_dEta[pc][1] = dN2_dEta(tableRow.xk[j]);
            dN_dEta[pc][2] = dN3_dEta(tableRow.xk[j]);
            dN_dEta[pc][3] = dN4_dEta(tableRow.xk[j]);
        }
    }
}

void UniversalElement::printShapeFunctionDerivatives()
{
    cout << "\ndNi/dKsi:" << endl;
    for (int pc = 0; pc < N * N; pc++) {
        for (int i = 0; i < 4; i++) {
            cout << dN_dKsi[pc][i] << "   " ;
        }
        cout << endl;
    }

    cout << "\ndNi/dEta:" << endl;
    for (int pc = 0; pc < N * N; pc++) {
        for (int i = 0; i < 4; i++) {
            cout << dN_dEta[pc][i] << "   ";
        }
        cout << endl;
    }
}

double dN1_dKsi(double eta) {
    return (-(1.0 / 4.0) * (1.0 - eta));
}

double dN2_dKsi(double eta) {
    return ((1.0 / 4.0) * (1.0 - eta));
}

double dN3_dKsi(double eta) {
    return ((1.0 / 4.0) * (1.0 + eta));
}

double dN4_dKsi(double eta) {
    return (-(1.0 / 4.0) * (1.0 + eta));
}

double dN1_dEta(double ksi) {
    return (-(1.0 / 4.0) * (1.0 - ksi));
}

double dN2_dEta(double ksi) {
    return (-(1.0 / 4.0) * (1.0 + ksi));
}

double dN3_dEta(double ksi) {
    return ((1.0 / 4.0) * (1.0 + ksi));
}

double dN4_dEta(double ksi) {
    return ((1.0 / 4.0) * (1.0 - ksi));
}