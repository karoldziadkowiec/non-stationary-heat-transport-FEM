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

void UniversalElement::calculateShapeFunctionDerivatives()
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int index = i * N + j;

            // dNi/dKsi
            dN_dKsi[index][0] = dN1_dKsi(tableRow.xk[i]);
            dN_dKsi[index][1] = dN2_dKsi(tableRow.xk[i]);
            dN_dKsi[index][2] = dN3_dKsi(tableRow.xk[i]);
            dN_dKsi[index][3] = dN4_dKsi(tableRow.xk[i]);

            // dNi/dEta
            dN_dEta[index][0] = dN1_dEta(tableRow.xk[j]);
            dN_dEta[index][1] = dN2_dEta(tableRow.xk[j]);
            dN_dEta[index][2] = dN3_dEta(tableRow.xk[j]);
            dN_dEta[index][3] = dN4_dEta(tableRow.xk[j]);
        }
    }
}

void UniversalElement::printShapeFunctionDerivatives()
{
    cout << "\ndNi/dKsi:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << dN_dKsi[i][j] << endl;
        }
    }

    cout << "\ndNi/dEta:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << dN_dEta[i][j] << endl;
        }
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