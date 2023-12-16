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

    for (int i = 0; i < 4; i++) {
        surface[i].ksiEtaMatrix = new double* [N];
        for (int j = 0; j < N; j++) {
            surface[i].ksiEtaMatrix[j] = new double[2] {};
        }
    }

    for (int i = 0; i < 4; i++) {
        surface[i].Ni = new double* [N];
        for (int j = 0; j < N; j++) {
            surface[i].Ni[j] = new double[4] {};
            for (int k = 0; k < 4; k++) {
                surface[i].Ni[j][k] = 0.0;
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        surface[i].Hbc_AtPci = new double** [N];
        for (int j = 0; j < N; j++) {
            surface[i].Hbc_AtPci[j] = new double* [4];
            for (int k = 0; k < 4; k++) {
                surface[i].Hbc_AtPci[j][k] = new double[4] {};
            }
        }
    }

    Ni_MatrixC = new double* [N * N];

    for (int i = 0; i < N * N; i++) {
        Ni_MatrixC[i] = new double[4] {};
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

/////////////////////////// H /////////////////////////////////////////////////////////////////////////////////////////////

void UniversalElement::calculateShapeFunctionDerivatives() // Pochodne funkcji kszta³tu dN/Kdsi oraz dN/dEta
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N); // Zwracanie wiersza zawieraj¹cego pc oraz wagi z tabeli kwadratury Gaussa w zale¿noœci od N

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

void UniversalElement::printShapeFunctionDerivatives() // Wypisanie pochodnych funkcji kszta³tu dN/Kdsi oraz dN/dEta
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

/////////////////////////// Hbc /////////////////////////////////////////////////////////////////////////////////////////////

void UniversalElement::calculateKsiEtaMatrix_Values() // Obliczenie Ksi oraz Eta dla funkcji kszat³tu  
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    //Surface 1
    for (int i = 0; i < N; i++) {
        surface[0].ksiEtaMatrix[i][0] = tableRow.xk[i];
        surface[0].ksiEtaMatrix[i][1] = -1;
    }

    //Surface 2
    for (int i = 0; i < N; i++) {
        surface[1].ksiEtaMatrix[i][0] = 1;
        surface[1].ksiEtaMatrix[i][1] = tableRow.xk[i];
    }

    //Surface 3
    for (int i = 0, j = N - 1; i < N; i++, j--) {
        surface[2].ksiEtaMatrix[i][0] = tableRow.xk[j];
        surface[2].ksiEtaMatrix[i][1] = 1;
    }

    //Surface 4
    for (int i = 0, j = N - 1; i < N; i++, j--) {
        surface[3].ksiEtaMatrix[i][0] = -1;
        surface[3].ksiEtaMatrix[i][1] = tableRow.xk[j];
    }
}

void UniversalElement::printKsiEtaMatrix_Values()
{
    for (int s = 0; s < 4; s++) {
        cout << "\nKsi and Eta values for surface " << s + 1 << endl;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < 2; j++) {
                cout << surface[s].ksiEtaMatrix[i][j] << "   ";
            }
            cout << endl;
        }
        cout << endl;
    }
}

void UniversalElement::calculateMatrixOfN_Values(int surf) // Obliczanie funkcji kszta³tu dla konkrtnej œciany
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0; i < N; i++) {
        surface[surf].Ni[i][0] = N_Function(surface[surf].ksiEtaMatrix[i][0], surface[surf].ksiEtaMatrix[i][1], 0);
        surface[surf].Ni[i][1] = N_Function(surface[surf].ksiEtaMatrix[i][0], surface[surf].ksiEtaMatrix[i][1], 1);
        surface[surf].Ni[i][2] = N_Function(surface[surf].ksiEtaMatrix[i][0], surface[surf].ksiEtaMatrix[i][1], 2);
        surface[surf].Ni[i][3] = N_Function(surface[surf].ksiEtaMatrix[i][0], surface[surf].ksiEtaMatrix[i][1], 3);
    }
}

/////////////////////////// C ///////////////////////////////////////////////////////////////////////////////////////////////

void UniversalElement::calculateMatrixOfN_ValuesMatrixC() // Obliczanie funkcji kszta³tu
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int pc = i * N + j;

            Ni_MatrixC[pc][0] = N_Function(tableRow.xk[j], tableRow.xk[i], 0);
            Ni_MatrixC[pc][1] = N_Function(tableRow.xk[j], tableRow.xk[i], 1);
            Ni_MatrixC[pc][2] = N_Function(tableRow.xk[j], tableRow.xk[i], 2);
            Ni_MatrixC[pc][3] = N_Function(tableRow.xk[j], tableRow.xk[i], 3);
        }
    }
}

void UniversalElement::printMatrixOfN_ValuesMatrixC()
{
    cout << "\nNi values for Matrix C:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << Ni_MatrixC[i][j] << "  ";
        }
        cout << endl;
    }
}

double N_Function(double ksi, double eta, int i) // Funkcje kszta³tu N1, N2, N3, N4
{
    switch (i)
    {
    case 0:
        return 0.25 * (1.0 - ksi) * (1.0 - eta);
    case 1:
        return 0.25 * (1.0 + ksi) * (1.0 - eta);
    case 2:
        return 0.25 * (1.0 + ksi) * (1.0 + eta);
    case 3:
        return 0.25 * (1.0 - ksi) * (1.0 + eta);
    default:
        cerr << "Invalid value of i in N function." << endl;
        return 0.0;
    }
}

void UniversalElement::printMatrixOfN_Values(int surf)
{
    cout << "\nSurface " << surf +1 << endl;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << surface[surf].Ni[i][j] << "\t";
        }
        cout << endl;
    }
    cout << endl;
}

// Pochodne funkcji kszta³tu
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