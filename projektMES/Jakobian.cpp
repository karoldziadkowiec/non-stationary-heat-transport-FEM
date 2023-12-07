#include <cmath>

#include "Grid.h"
#include "Jakobian.h"

void Jakobian::calculateDerivativesAtPci(const UniversalElement& universalElement, const Grid& grid, int elementNumber, int pc)
{
    dx_dKsi = universalElement.dN_dKsi[pc][0] * grid.nodes[grid.elements[elementNumber].id[2] - 1].x +
        universalElement.dN_dKsi[pc][1] * grid.nodes[grid.elements[elementNumber].id[3] - 1].x +
        universalElement.dN_dKsi[pc][2] * grid.nodes[grid.elements[elementNumber].id[0] - 1].x +
        universalElement.dN_dKsi[pc][3] * grid.nodes[grid.elements[elementNumber].id[1] - 1].x;

    dy_dKsi = universalElement.dN_dKsi[pc][0] * grid.nodes[grid.elements[elementNumber].id[2] - 1].y +
        universalElement.dN_dKsi[pc][1] * grid.nodes[grid.elements[elementNumber].id[3] - 1].y +
        universalElement.dN_dKsi[pc][2] * grid.nodes[grid.elements[elementNumber].id[0] - 1].y +
        universalElement.dN_dKsi[pc][3] * grid.nodes[grid.elements[elementNumber].id[1] - 1].y;

    dx_dEta = universalElement.dN_dEta[pc][0] * grid.nodes[grid.elements[elementNumber].id[2] - 1].x +
        universalElement.dN_dEta[pc][1] * grid.nodes[grid.elements[elementNumber].id[3] - 1].x +
        universalElement.dN_dEta[pc][2] * grid.nodes[grid.elements[elementNumber].id[0] - 1].x +
        universalElement.dN_dEta[pc][3] * grid.nodes[grid.elements[elementNumber].id[1] - 1].x;

    dy_dEta = universalElement.dN_dEta[pc][0] * grid.nodes[grid.elements[elementNumber].id[2] - 1].y +
        universalElement.dN_dEta[pc][1] * grid.nodes[grid.elements[elementNumber].id[3] - 1].y +
        universalElement.dN_dEta[pc][2] * grid.nodes[grid.elements[elementNumber].id[0] - 1].y +
        universalElement.dN_dEta[pc][3] * grid.nodes[grid.elements[elementNumber].id[1] - 1].y;
}

void Jakobian::printJakobianMatrix()
{
    cout << "\nJakobian Matrix" << endl;
    cout << dx_dKsi << "\t" << dy_dKsi << endl;
    cout << dx_dEta << "\t" << dy_dEta << endl;
}

double Jakobian::calculateDetJ()
{
    double detJ = (dx_dKsi * dy_dEta) - (dy_dKsi * dx_dEta);
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

    Hbc_AtPci = new double** [4];
    Hbci = new double** [4];

    Cpci = new double** [N * N];

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
        Cpci[i] = new double* [N * N];

        for (int j = 0; j < 4; j++) {
            Hpci[i][j] = new double[4] {};
            Cpci[i][j] = new double[4] {};
        }
    }

    for (int i = 0; i < 4; i++) {
        Hbc_AtPci[i] = new double* [4];
        for (int j = 0; j < 4; j++) {
            Hbc_AtPci[i][j] = new double[4] {};
        }
    }

    for (int i = 0; i < 4; i++) {
        Hbci[i] = new double* [4];
        for (int j = 0; j < 4; j++) {
            Hbci[i][j] = new double[4] {};
        }
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
    }
    delete[] dN_dx;
    delete[] dN_dy;

    delete[] xH_AtPci;
    delete[] yH_AtPci;

    delete[] Hpci;
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

void Jakobian::calculateMatrixHForXandYForPci(int pc)
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

void Jakobian::calculateMatrixHpci(int pc, int conductivity)
{
    int kt = conductivity; // conductivity
    dV = calculate1_DetJ(); //volume of the integrated element 

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Hpci[pc][i][j] = kt * (xH_AtPci[pc][i][j] + yH_AtPci[pc][i][j]) * dV;
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

void Jakobian::calculateMatrixH(const Grid& grid, int elementNumber)
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            grid.elements[elementNumber].H[i][j] = 0;

            for (int pc = 0; pc < N * N; pc++)
            {
                int w1 = pc / N;
                int w2 = pc % N;

                grid.elements[elementNumber].H[i][j] += Hpci[pc][i][j] * tableRow.wk[w1] * tableRow.wk[w2];
            }
        }
    }
}

void Jakobian::printMatrixH(const Grid& grid, int elementNumber)
{
    cout << "\nMatrix [H]:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << grid.elements[elementNumber].H[i][j] << "   ";
        }
        cout << endl;
    }
}
/////////////////////////////////// Hbc //////////////////////////////////////////////
double Jakobian::calculateHbcDetJ(const Grid& grid, int Nx, int Nk)
{
    double L = sqrt(pow((grid.nodes[Nx].x - grid.nodes[Nk].x), 2) + pow((grid.nodes[Nx].y - grid.nodes[Nk].y), 2));
    double DetJ = L / 2;
    return DetJ;
}

void Jakobian::printHbcDetJ(const Grid& grid, int Nx, int Nk)
{
    cout << "DetJ for surface = " << calculateHbcDetJ(grid, Nx, Nk) << endl;
}

void Jakobian::calculateMatrixHbciForPci(const UniversalElement& universalElement, int surface)
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0, n = N - 1; i < N; i++, n--)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                if (surface == 0 || surface == 1) {
                    Hbc_AtPci[surface][j][k] += tableRow.wk[i] * universalElement.surface[surface].Ni[i][j] * universalElement.surface[surface].Ni[i][k];
                }
                else {
                    Hbc_AtPci[surface][j][k] += tableRow.wk[n] * universalElement.surface[surface].Ni[i][j] * universalElement.surface[surface].Ni[i][k];
                }
            }
        }
    }
}

void Jakobian::calculateMatrixHbci(int surface, int alfa, const Grid& grid, int Nx, int Nk)
{
    int a = alfa; // alfa
    dS = calculateHbcDetJ(grid, Nx, Nk); // area of the integrated element 

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Hbci[surface][i][j] = a * Hbc_AtPci[surface][i][j] * dS;
        }
    }
}

void Jakobian::printMatrixHbci()
{
    for (int s = 0; s < 4; s++) {
        cout << "\nHbc" << s + 1 << ":" << endl;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << Hbci[s][i][j] << "  ";
            }
            cout << endl;
        }
    }
}

void Jakobian::calculateMatrixHbc(const Grid& grid, int elementNumber)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            grid.elements[elementNumber].Hbc[i][j] = Hbci[0][i][j] + Hbci[1][i][j] + Hbci[2][i][j] + Hbci[3][i][j];
        }
    }
}

void Jakobian::printMatrixHbc(const Grid& grid, int elementNumber)
{
    cout << "\nMatrix [Hbc]:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << grid.elements[elementNumber].Hbc[i][j] << "              ";
        }
        cout << endl;
    }
}

void Jakobian::zeroMatrixHbci()
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            for (int k = 0; k < 4; k++)
            {
                Hbc_AtPci[i][j][k] = 0.0;
                Hbci[i][j][k] = 0.0;
            }
        }
    }
}

void Jakobian::sumMatrixH_Hbc(const Grid& grid, int elementNumber)
{
    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            grid.elements[elementNumber].H[i][j] += grid.elements[elementNumber].Hbc[i][j];
        }
    }
}

///////////////////////// Vector P ////////////////////////////////
void Jakobian::zeroVectorP(const Grid& grid, int elementNumber)
{
    for (int i = 0; i < 4; i++)
    {
        grid.elements[elementNumber].P[i] = 0.0;
    }
}

void Jakobian::calculateVectorP_ForPci(const Grid& grid, const UniversalElement& universalElement, int surface, int elementNumber, int tot, int alfa, int Nx, int Nk)
{
    int a = alfa; // alfa
    dS = calculateHbcDetJ(grid, Nx, Nk); // area of the integrated element
    int t = tot; // ambient temperature
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0, n = N - 1; i < N; i++, n--)
    {
        for (int j = 0; j < 4; j++)
        {
            if (surface == 0 || surface == 1) {
                grid.elements[elementNumber].P[j] += a * (tableRow.wk[i] * universalElement.surface[surface].Ni[i][j] * t) * dS;
            }
            else {
                grid.elements[elementNumber].P[j] += a * (tableRow.wk[n] * universalElement.surface[surface].Ni[i][j] * t) * dS;
            }
        }
    }
}

void Jakobian::printVectorP(const Grid& grid, int elementNumber)
{
    cout << "\nVector {P}:" << endl;
    for (int i = 0; i < 4; i++) {
        cout << grid.elements[elementNumber].P[i] << "   ";
    }
    cout << endl;
}

// Matrix C
void Jakobian::calculateMatrixCpci(const UniversalElement& universalElement, int pc, int specificHeat, int density)
{
    int c = specificHeat; // specific heat
    int ro = density; // density
    dV = calculate1_DetJ(); //volume of the integrated element 

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            Cpci[pc][i][j] = c * ro * (universalElement.Ni_MatrixC[pc][i] * universalElement.Ni_MatrixC[pc][j]) * dV;
        }
    }
}

void Jakobian::printMatrixCpci()
{
    for (int pc = 0; pc < N * N; pc++) {
        cout << "\nCpc" << pc + 1 << ":" << endl;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                cout << Cpci[pc][i][j] << "  ";
            }
            cout << endl;
        }
    }
}

void Jakobian::calculateMatrixC(const Grid& grid, int elementNumber)
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N);

    for (int i = 0; i < 4; i++)
    {
        for (int j = 0; j < 4; j++)
        {
            grid.elements[elementNumber].C[i][j] = 0;

            for (int pc = 0; pc < N * N; pc++)
            {
                int w1 = pc / N;
                int w2 = pc % N;

                grid.elements[elementNumber].C[i][j] += Cpci[pc][i][j] * tableRow.wk[w1] * tableRow.wk[w2];
            }
        }
    }
}

void Jakobian::printMatrixC(const Grid& grid, int elementNumber)
{
    cout << "\nMatrix [C]:" << endl;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            cout << grid.elements[elementNumber].C[i][j] << "   ";
        }
        cout << endl;
    }
}