#include <cmath>

#include "Grid.h"
#include "Jakobian.h"

void Jakobian::calculateDerivativesAtPci(const UniversalElement& universalElement, const Grid& grid, int i)
{
    dx_dKsi = universalElement.dN_dKsi[0][0] * grid.nodes[grid.elements[i].id[0] - 1].x +
        universalElement.dN_dKsi[0][1] * grid.nodes[grid.elements[i].id[1] - 1].x +
        universalElement.dN_dKsi[0][2] * grid.nodes[grid.elements[i].id[2] - 1].x +
        universalElement.dN_dKsi[0][3] * grid.nodes[grid.elements[i].id[3] - 1].x;

    dy_dKsi = -1 * (universalElement.dN_dKsi[0][0] * grid.nodes[grid.elements[i].id[0] - 1].y +
        universalElement.dN_dKsi[0][1] * grid.nodes[grid.elements[i].id[1] - 1].y +
        universalElement.dN_dKsi[0][2] * grid.nodes[grid.elements[i].id[2] - 1].y +
        universalElement.dN_dKsi[0][3] * grid.nodes[grid.elements[i].id[3] - 1].y);

    dx_dEta = -1 * (universalElement.dN_dEta[0][0] * grid.nodes[grid.elements[i].id[0] - 1].x +
        universalElement.dN_dEta[0][1] * grid.nodes[grid.elements[i].id[1] - 1].x +
        universalElement.dN_dEta[0][2] * grid.nodes[grid.elements[i].id[2] - 1].x +
        universalElement.dN_dEta[0][3] * grid.nodes[grid.elements[i].id[3] - 1].x);

    dy_dEta = universalElement.dN_dEta[0][0] * grid.nodes[grid.elements[i].id[0] - 1].y +
        universalElement.dN_dEta[0][1] * grid.nodes[grid.elements[i].id[1] - 1].y +
        universalElement.dN_dEta[0][2] * grid.nodes[grid.elements[i].id[2] - 1].y +
        universalElement.dN_dEta[0][3] * grid.nodes[grid.elements[i].id[3] - 1].y;
}

void Jakobian::printJakobianMatrix()
{
    cout << "Jakobian Matrix" << endl;
    cout << dy_dEta << " " << dy_dKsi << endl;
    cout << dx_dEta << " " << dx_dKsi << endl;
}

double Jakobian::calculateDetJ()
{
    double detJ = (dy_dEta * dx_dKsi) - (dy_dKsi * dx_dEta);
    return detJ;
}

double Jakobian::calculate1_DetJ()
{
    double result = 1.0 / calculateDetJ();
    return result;
}

void Jakobian::calculateJakobianMatrix()
{
    double value1_DetJ = calculate1_DetJ();
    dx_dKsi = value1_DetJ * dx_dKsi;
    dy_dKsi = value1_DetJ * dy_dKsi;
    dx_dEta = value1_DetJ * dx_dEta;
    dy_dEta = value1_DetJ * dy_dEta;
}

Jakobian::Jakobian(int N)
{
    this->N = N;

    dN_dx = new double* [N * N];
    dN_dy = new double* [N * N];

    xH_AtPc1 = new double* [N * N];
    yH_AtPc1 = new double* [N * N];
    xH_AtPc2 = new double* [N * N];
    yH_AtPc2 = new double* [N * N];
    xH_AtPc3 = new double* [N * N];
    yH_AtPc3 = new double* [N * N];
    xH_AtPc4 = new double* [N * N];
    yH_AtPc4 = new double* [N * N];

    Hpc1 = new double* [N * N];
    Hpc2 = new double* [N * N];
    Hpc3 = new double* [N * N];
    Hpc4 = new double* [N * N];

    H = new double* [N * N];

    for (int i = 0; i < N * N; i++) {
        dN_dx[i] = new double[4] {};
        dN_dy[i] = new double[4] {};

        xH_AtPc1[i] = new double[4] {};
        yH_AtPc1[i] = new double[4] {};
        xH_AtPc2[i] = new double[4] {};
        yH_AtPc2[i] = new double[4] {};
        xH_AtPc3[i] = new double[4] {};
        yH_AtPc3[i] = new double[4] {};
        xH_AtPc4[i] = new double[4] {};
        yH_AtPc4[i] = new double[4] {};

        Hpc1[i] = new double[4] {};
        Hpc2[i] = new double[4] {};
        Hpc3[i] = new double[4] {};
        Hpc4[i] = new double[4] {};

        H[i] = new double[4] {};
    }
}

Jakobian::~Jakobian()
{
    for (int i = 0; i < N * N; i++) {
        delete[] dN_dx[i];
        delete[] dN_dy[i];

        delete[] xH_AtPc1[i];
        delete[] yH_AtPc1[i];
        delete[] xH_AtPc2[i];
        delete[] yH_AtPc2[i];
        delete[] xH_AtPc3[i];
        delete[] yH_AtPc3[i];
        delete[] xH_AtPc4[i];
        delete[] yH_AtPc4[i];

        delete[] Hpc1[i];
        delete[] Hpc2[i];
        delete[] Hpc3[i];
        delete[] Hpc4[i];

        delete[] H[i];
    }
    delete[] dN_dx;
    delete[] dN_dy;

    delete[] xH_AtPc1;
    delete[] yH_AtPc1;
    delete[] xH_AtPc2;
    delete[] yH_AtPc2;
    delete[] xH_AtPc3;
    delete[] yH_AtPc3;
    delete[] xH_AtPc4;
    delete[] yH_AtPc4;

    delete[] Hpc1;
    delete[] Hpc2;
    delete[] Hpc3;
    delete[] Hpc4;

    delete[] H;
}

void Jakobian::calculateShapeFunctionDerivativesForPci(const UniversalElement& universalElement)
{
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int index = i * N + j;

             dN_dx[index][0] = dy_dEta * universalElement.dN_dKsi[index][0] + dy_dKsi * universalElement.dN_dEta[index][0];
             dN_dx[index][1] = dy_dEta * universalElement.dN_dKsi[index][1] + dy_dKsi * universalElement.dN_dEta[index][1];
             dN_dx[index][2] = dy_dEta * universalElement.dN_dKsi[index][2] + dy_dKsi * universalElement.dN_dEta[index][2];
             dN_dx[index][3] = dy_dEta * universalElement.dN_dKsi[index][3] + dy_dKsi * universalElement.dN_dEta[index][3];

             dN_dy[index][0] = dx_dEta * universalElement.dN_dKsi[index][0] + dx_dKsi * universalElement.dN_dEta[index][0];
             dN_dy[index][1] = dx_dEta * universalElement.dN_dKsi[index][1] + dx_dKsi * universalElement.dN_dEta[index][1];
             dN_dy[index][2] = dx_dEta * universalElement.dN_dKsi[index][2] + dx_dKsi * universalElement.dN_dEta[index][2];
             dN_dy[index][3] = dx_dEta * universalElement.dN_dKsi[index][3] + dx_dKsi * universalElement.dN_dEta[index][3];
        }
    }
}


void Jakobian::printShapeFunctionDerivativesForPci()
{
    cout << "\ndNi/dx:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << dN_dx[i][j] << endl;
        }
    }

    cout << "\ndNi/dy:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << dN_dy[i][j] << endl;
        }
    }
}

void Jakobian::calculateMatrixHForXandYForPci()
{
    for (int i = 0; i < N * N; ++i)
    {
        for (int j = 0; j < N * N; ++j)
        {
            xH_AtPc1[i][j] = dN_dx[0][i] * dN_dx[0][j];
            yH_AtPc1[i][j] = dN_dy[0][i] * dN_dy[0][j];

            xH_AtPc2[i][j] = dN_dx[1][i] * dN_dx[1][j];
            yH_AtPc2[i][j] = dN_dy[1][i] * dN_dy[1][j];

            xH_AtPc3[i][j] = dN_dx[2][i] * dN_dx[2][j];
            yH_AtPc3[i][j] = dN_dy[2][i] * dN_dy[2][j];

            xH_AtPc4[i][j] = dN_dx[3][i] * dN_dx[3][j];
            yH_AtPc4[i][j] = dN_dy[3][i] * dN_dy[3][j];
        }
    }
}

void Jakobian::printMatrixHForXandYForPci()
{
    cout << "\nMatrix H in pc1 for x:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << xH_AtPc1[i][j] << endl;
        }
    }

    cout << "\nMatrix H in pc1 for y:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << yH_AtPc1[i][j] << endl;
        }
    }

    cout << "\nMatrix H in pc2 for x:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << xH_AtPc2[i][j] << endl;
        }
    }

    cout << "\nMatrix H in pc2 for y:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << yH_AtPc2[i][j] << endl;
        }
    }

    cout << "\nMatrix H in pc3 for x:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << xH_AtPc3[i][j] << endl;
        }
    }

    cout << "\nMatrix H in pc3 for y:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << yH_AtPc3[i][j] << endl;
        }
    }

    cout << "\nMatrix H in pc4 for x:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << xH_AtPc4[i][j] << endl;
        }
    }

    cout << "\nMatrix H in pc4 for y:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << yH_AtPc4[i][j] << endl;
        }
    }
}

void Jakobian::calculateMatrixHpci(int kt)
{
    dV = calculate1_DetJ(); //area of the integrated element 

    for (int i = 0; i < N * N; ++i)
    {
        for (int j = 0; j < N * N; ++j)
        {
            Hpc1[i][j] = kt * (xH_AtPc1[i][j] + yH_AtPc1[i][j]) * dV;

            Hpc2[i][j] = kt * (xH_AtPc2[i][j] + yH_AtPc2[i][j]) * dV;

            Hpc3[i][j] = kt * (xH_AtPc3[i][j] + yH_AtPc3[i][j]) * dV;

            Hpc4[i][j] = kt * (xH_AtPc4[i][j] + yH_AtPc4[i][j]) * dV;
        }
    }
}

void Jakobian::printMatrixHpci()
{
    cout << "\nHpc1:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << Hpc1[i][j] << endl;
        }
    }

    cout << "\nHpc2:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << Hpc2[i][j] << endl;
        }
    }

    cout << "\nHpc3:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << Hpc3[i][j] << endl;
        }
    }

    cout << "\nHpc4:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << Hpc4[i][j] << endl;
        }
    }
}

void Jakobian::calculateMatrixH()
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    double w1 = tableRow.wk[0];
    double w2 = tableRow.wk[1];

    for (int i = 0; i < N * N; ++i)
    {
        for (int j = 0; j < N * N; ++j)
        {
            H[i][j] = (Hpc1[i][j] * w1 * w1) +
                (Hpc2[i][j] * w2 * w1) +
                (Hpc3[i][j] * w1 * w2) +
                (Hpc4[i][j] * w2 * w2);
        }
    }
}

void Jakobian::printMatrixH()
{
    cout << "\nMatrix [H]:" << endl;
    for (int i = 0; i < N * N; i++) {
        for (int j = 0; j < 4; j++) {
            cout << "[" << i + 1 << "][" << j + 1 << "] = " << H[i][j] << endl;
        }
    }
}