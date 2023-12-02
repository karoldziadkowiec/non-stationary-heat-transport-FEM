#ifndef AGREGATION_H
#define AGREGATION_H
#include <iostream>
#include <cmath>

#include "GaussQuadrature.h"
#include "UniversalElement.h"

using namespace std;

// System of Equation
struct SoE {
	int n;
	double** HG;
	double** CG;
	double* P;
	double* t;

	SoE(int n);
	~SoE();

	void aggregateMatrixH(const Grid& grid, int elementNumber);
	void printAggregatedMatrixH();
	void aggregateVectorP(const Grid& grid, int elementNumber);
	void printAggregatedVectorP();
	void solveSoE();
	void printSoE();
};

#endif