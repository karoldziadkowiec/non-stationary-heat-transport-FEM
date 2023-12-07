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
	double* PG;
	double* t;

	SoE(int n);
	~SoE();

	void aggregateMatrixHG(const Grid& grid, int elementNumber);
	void printAggregatedMatrixHG();

	void aggregateVectorPG(const Grid& grid, int elementNumber);
	void printAggregatedVectorPG();

	void aggregateMatrixCG(const Grid& grid, int elementNumber);
	void printAggregatedMatrixCG();

	void calculateMatrixHplusC_dT(const Grid& grid, int elementsNumber, int dt);
	void printMatrixHplusC_dT(const Grid& grid, int elementsNumber);
	void zeroAggregatedMatrixH();
	void aggregateMatrixHplusC_dT(const Grid& grid, int elementsNumber);
	void calculateMatrixCt0_dTplusP(const Grid& grid, int elementsNumber, int dt, int t);
	void printMatrixCt0_dTplusP(const Grid& grid, int elementsNumber);
	void zeroAggregatedMatrixP();
	void aggregateMatrixCt0_dTplusP(const Grid& grid, int elementsNumber);
	void displayMinMaxTemperature(int dt);

	void solveSoE();
	void printSoE();
};

#endif