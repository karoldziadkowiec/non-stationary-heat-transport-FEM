#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <cmath>

#include "GaussQuadrature.h"

GaussQuadrature returnRowOfGaussTable(int N) // Tabela kwadratury Gaussa zawieraj¹ca info o pc oraz wagach
{
    GaussQuadrature tableRow;

    switch (N)
    {
    case 2:
        tableRow.k = new int[2] {0, 1};
        tableRow.xk = new double[2] {-1.0/sqrt(3), 1.0/sqrt(3)};
        tableRow.wk = new double[2] {1.0, 1.0};
        break;
    case 3:
        tableRow.k = new int[3] {0, 1, 2};
        tableRow.xk = new double[3] {-sqrt(3.0/5.0), 0.0, sqrt(3.0/5.0)};
        tableRow.wk = new double[3] {5.0/9.0, 8.0/9.0, 5.0/9.0};
        break;
    case 4:
        tableRow.k = new int[4] {0, 1, 2, 3};
        tableRow.xk = new double[4] {-0.861136, -0.339981, 0.339981, 0.861136};
        tableRow.wk = new double[4] {0.347855, 0.652145, 0.652145, 0.347855};
        break;
    default:
        cerr << "Unknown number of nodes." << endl;
        tableRow.k = nullptr;
        tableRow.xk = nullptr;
        tableRow.wk = nullptr;
    }
    return tableRow;
}

double gauss(int N, double (*function)(double))
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N); //Returning weights and Gaussian points for given number of points N

    double integral = 0.0; // Integral result

    for (int i = 0; i < N; i++) {
        integral += tableRow.wk[i] * function(tableRow.xk[i]);
    }
    return integral;
}

double gauss(int N, double (*function)(double, double))
{
    GaussQuadrature tableRow = returnRowOfGaussTable(N); //Returning weights and Gaussian points for given number of points N

    double integral = 0.0; // Integral result

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            integral += tableRow.wk[i] * tableRow.wk[j] * function(tableRow.xk[i], tableRow.xk[j]);
        }
    }
    return integral;
}