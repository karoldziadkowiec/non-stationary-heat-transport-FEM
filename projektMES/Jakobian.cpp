#include <cmath>

#include "Grid.h"
#include "Jakobian.h"

void Jakobian::calculateDerivativesAtPci(const UniversalElement& universalElement, const Grid& grid, int elementNumber, int pc)
{
    dx_dKsi = universalElement.dN_dKsi[pc][0] * grid.nodes[grid.elements[elementNumber].id[0] - 1].x +
        universalElement.dN_dKsi[pc][1] * grid.nodes[grid.elements[elementNumber].id[1] - 1].x +
        universalElement.dN_dKsi[pc][2] * grid.nodes[grid.elements[elementNumber].id[2] - 1].x +
        universalElement.dN_dKsi[pc][3] * grid.nodes[grid.elements[elementNumber].id[3] - 1].x;

    dy_dKsi = universalElement.dN_dKsi[pc][0] * grid.nodes[grid.elements[elementNumber].id[0] - 1].y +
        universalElement.dN_dKsi[pc][1] * grid.nodes[grid.elements[elementNumber].id[1] - 1].y +
        universalElement.dN_dKsi[pc][2] * grid.nodes[grid.elements[elementNumber].id[2] - 1].y +
        universalElement.dN_dKsi[pc][3] * grid.nodes[grid.elements[elementNumber].id[3] - 1].y;

    dx_dEta = universalElement.dN_dEta[pc][0] * grid.nodes[grid.elements[elementNumber].id[0] - 1].x +
        universalElement.dN_dEta[pc][1] * grid.nodes[grid.elements[elementNumber].id[1] - 1].x +
        universalElement.dN_dEta[pc][2] * grid.nodes[grid.elements[elementNumber].id[2] - 1].x +
        universalElement.dN_dEta[pc][3] * grid.nodes[grid.elements[elementNumber].id[3] - 1].x;

    dy_dEta = universalElement.dN_dEta[pc][0] * grid.nodes[grid.elements[elementNumber].id[0] - 1].y +
        universalElement.dN_dEta[pc][1] * grid.nodes[grid.elements[elementNumber].id[1] - 1].y +
        universalElement.dN_dEta[pc][2] * grid.nodes[grid.elements[elementNumber].id[2] - 1].y +
        universalElement.dN_dEta[pc][3] * grid.nodes[grid.elements[elementNumber].id[3] - 1].y;
}

void Jakobian::printJakobianMatrix()
{
    cout << "\nJakobian Matrix" << endl;
    cout << dx_dKsi << "\t" << dy_dKsi << endl;
    cout << dx_dEta << "\t" << dy_dEta << endl;
}

double Jakobian::calculateDetJ()
{
    double detJ = (dy_dEta * dx_dKsi) - (dy_dKsi * dx_dEta);
    return detJ;
}

void Jakobian::printDetJ()
{
    cout << "\nDetJ = " << calculateDetJ() << endl;
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

    xH_AtPci = new double** [N * N];
    yH_AtPci = new double** [N * N];

    Hpci = new double** [N * N];

    H = new double* [N * N];

    for (int i = 0; i < N * N; i++) {
        dN_dx[i] = new double[4] {};
        dN_dy[i] = new double[4] {};

        xH_AtPci[i] = new double* [N * N];
        yH_AtPci[i] = new double* [N * N];

        for (int j = 0; j < 4; j++) {
            xH_AtPci[i][j] = new double[4] {};
            yH_AtPci[i][j] = new double[4] {};
        }

        Hpci[i] = new double* [N * N];

        for (int j = 0; j < 4; j++) {
            Hpci[i][j] = new double[4] {};
        }

        H[i] = new double[4] {};
    }
}

Jakobian::~Jakobian()
{
    for (int i = 0; i < N * N; i++) {
        delete[] dN_dx[i];
        delete[] dN_dy[i];

        for (int j = 0; j < 4; j++) {
            delete[] xH_AtPci[i][j];
            delete[] yH_AtPci[i][j];
        }

        delete[] xH_AtPci[i];
        delete[] yH_AtPci[i];

        delete[] Hpci[i];

        delete[] H[i];
    }
    delete[] dN_dx;
    delete[] dN_dy;

    delete[] xH_AtPci;
    delete[] yH_AtPci;

    delete[] Hpci;

    delete[] H;
}

void Jakobian::calculateShapeFunctionDerivativesForPci(const UniversalElement& universalElement, int pc)
{
    dN_dx[pc][0] = dy_dEta * universalElement.dN_dKsi[pc][0] + (-1 * dy_dKsi) * universalElement.dN_dEta[pc][0];
    dN_dx[pc][1] = dy_dEta * universalElement.dN_dKsi[pc][1] + (-1 * dy_dKsi) * universalElement.dN_dEta[pc][1];
    dN_dx[pc][2] = dy_dEta * universalElement.dN_dKsi[pc][2] + (-1 * dy_dKsi) * universalElement.dN_dEta[pc][2];
    dN_dx[pc][3] = dy_dEta * universalElement.dN_dKsi[pc][3] + (-1 * dy_dKsi) * universalElement.dN_dEta[pc][3];

    dN_dy[pc][0] = (-1 * dx_dEta) * universalElement.dN_dKsi[pc][0] + dx_dKsi * universalElement.dN_dEta[pc][0];
    dN_dy[pc][1] = (-1 * dx_dEta) * universalElement.dN_dKsi[pc][1] + dx_dKsi * universalElement.dN_dEta[pc][1];
    dN_dy[pc][2] = (-1 * dx_dEta) * universalElement.dN_dKsi[pc][2] + dx_dKsi * universalElement.dN_dEta[pc][2];
    dN_dy[pc][3] = (-1 * dx_dEta) * universalElement.dN_dKsi[pc][3] + dx_dKsi * universalElement.dN_dEta[pc][3];
}


void Jakobian::printShapeFunctionDerivativesForPci(int pc)
{
    cout << "\ndNi/dx:" << endl;
    for (int i = 0; i < 4; i++) {
        cout << dN_dx[pc][i] << "  ";
    }
    cout << endl;

    cout << "\ndNi/dy:" << endl;
    for (int i = 0; i < 4; i++) {
        cout << dN_dy[pc][i] << "  ";
    }
    cout << endl;
}

void Jakobian::calculateMatrixHForXandYForPci()
{
    for (int pc = 0; pc < N * N; pc++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                xH_AtPci[pc][i][j] = dN_dx[pc][i] * dN_dx[pc][j];
                yH_AtPci[pc][i][j] = dN_dy[pc][i] * dN_dy[pc][j];
            }
        }
    }
}

void Jakobian::printMatrixHForXandYForPci()
{
    for (int pc = 0; pc < N * N; pc++) {
        cout << "\nMatrix H in pc" << pc + 1 << " for dx:" << endl;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << xH_AtPci[pc][i][j] << "  ";
            }
            cout << endl;
        }

        cout << "\nMatrix H in pc" << pc + 1 << " for dy:" << endl;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << yH_AtPci[pc][i][j] << "  ";
            }
            cout << endl;
        }
    }
}

void Jakobian::calculateMatrixHpci(int conductivity)
{
    int kt = conductivity; // conductivity
    dV = calculate1_DetJ(); //area of the integrated element 

    for (int pc = 0; pc < N * N; pc++)
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                Hpci[pc][i][j] = kt * (xH_AtPci[pc][i][j] + yH_AtPci[pc][i][j]) * dV;
            }
        }
    }
}

void Jakobian::printMatrixHpci()
{
    for (int pc = 0; pc < N * N; pc++) {
        cout << "\nHpc" << pc + 1 << ":" << endl;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << Hpci[pc][i][j] << "  ";
            }
            cout << endl;
        }
    }
}

void Jakobian::calculateMatrixH()
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            H[i][j] = 0;

            for (int pc = 0; pc < N * N; pc++)
            {
                int w1 = pc / N;
                int w2 = pc % N;

                H[i][j] += Hpci[pc][i][j] * tableRow.wk[w1] * tableRow.wk[w2];
            }
        }
    }
}

void Jakobian::printMatrixH()
{
    cout << "\nMatrix [H]:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << H[i][j] << "   ";
        }
        cout << endl;
    }
}