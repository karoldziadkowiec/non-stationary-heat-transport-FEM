#include <cmath>

#include "Grid.h"
#include "Agregation.h"

SoE::SoE(int n)
{
	this->n = n;

	HG = new double* [n];
	CG = new double* [n];
	for (int i = 0; i < n; i++) {
		HG[i] = new double[n] {};
		CG[i] = new double[n] {};
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			HG[i][j] = 0.0;
			CG[i][j] = 0.0;
		}
	}

	PG = new double[n];
	t = new double[n];
	for (int i = 0; i < n; i++) {
		PG[i] = 0.0;
		t[i] = 0.0;
	}
}

SoE::~SoE()
{
	for (int i = 0; i < n; i++) {
		delete[] HG[i];
		delete[] CG[i];
	}
	delete[] HG;
	delete[] CG;
	delete[] PG;
}

void SoE::aggregateMatrixHG(const Grid& grid, int elementNumber)
{
	HG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].H[0][0];
	HG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].H[0][1];
	HG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].H[0][2];
	HG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].H[0][3];

	HG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].H[1][0];
	HG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].H[1][1];
	HG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].H[1][2];
	HG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].H[1][3];

	HG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].H[2][0];
	HG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].H[2][1];
	HG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].H[2][2];
	HG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].H[2][3];

	HG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].H[3][0];
	HG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].H[3][1];
	HG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].H[3][2];
	HG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].H[3][3];
}

void SoE::printAggregatedMatrixHG()
{
	cout << "\nAgrregated Matrix [HG] n x n:" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << HG[i][j] << " ";
		}
		cout << endl;
	}
}

void SoE::aggregateVectorPG(const Grid& grid, int elementNumber)
{
	PG[grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].P[0];
	PG[grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].P[1];
	PG[grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].P[2];
	PG[grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].P[3];
}

void SoE::printAggregatedVectorPG()
{
	cout << "\nAgrregated Vector {P} n x 1:" << endl;
	for (int i = 0; i < n; i++) {
		cout << PG[i] << "  ";
	}
	cout << endl;
}

void SoE::aggregateMatrixCG(const Grid& grid, int elementNumber)
{
	CG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].C[0][0];
	CG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].C[0][1];
	CG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].C[0][2];
	CG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].C[0][3];

	CG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].C[1][0];
	CG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].C[1][1];
	CG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].C[1][2];
	CG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].C[1][3];

	CG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].C[2][0];
	CG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].C[2][1];
	CG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].C[2][2];
	CG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].C[2][3];

	CG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].C[3][0];
	CG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].C[3][1];
	CG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].C[3][2];
	CG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].C[3][3];
}

void SoE::printAggregatedMatrixCG()
{
	cout << "\nAgrregated Matrix [CG] n x n:" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << CG[i][j] << " ";
		}
		cout << endl;
	}
}

void SoE::calculateMatrixHplusC_dT(const Grid& grid, int elementsNumber, int dt)
{
	int deltaTau = dt;
	for (int elNumber = 0; elNumber < elementsNumber; elNumber++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				grid.elements[elNumber].HplusC_dT[i][j] = grid.elements[elNumber].H[i][j] + (grid.elements[elNumber].C[i][j]/deltaTau);
			}
		}
	}
}

void SoE::printMatrixHplusC_dT(const Grid& grid, int elementsNumber)
{
	for (int elNumber = 0; elNumber < elementsNumber; elNumber++) {
		cout << "\nMatrix " << elNumber + 1 << " [H + C/dT]:" << endl;
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				cout << grid.elements[elNumber].HplusC_dT[i][j] << " ";
			}
			cout << endl;
		}
	}
}

void SoE::zeroAggregatedMatrixH()
{
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			HG[i][j] = 0.0;
		}
	}
}

void SoE::aggregateMatrixHplusC_dT(const Grid& grid, int elementsNumber)
{
	for (int elementNumber = 0; elementNumber < elementsNumber; elementNumber++) {
		HG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].HplusC_dT[0][0];
		HG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].HplusC_dT[0][1];
		HG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].HplusC_dT[0][2];
		HG[grid.elements[elementNumber].id[2] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].HplusC_dT[0][3];

		HG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].HplusC_dT[1][0];
		HG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].HplusC_dT[1][1];
		HG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].HplusC_dT[1][2];
		HG[grid.elements[elementNumber].id[3] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].HplusC_dT[1][3];

		HG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].HplusC_dT[2][0];
		HG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].HplusC_dT[2][1];
		HG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].HplusC_dT[2][2];
		HG[grid.elements[elementNumber].id[0] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].HplusC_dT[2][3];

		HG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].HplusC_dT[3][0];
		HG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].HplusC_dT[3][1];
		HG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].HplusC_dT[3][2];
		HG[grid.elements[elementNumber].id[1] - 1][grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].HplusC_dT[3][3];
	}
}

void SoE::calculateMatrixCt0_dTplusP(const Grid& grid, int elementsNumber, int dt)
{
	int deltaTau = dt;

	for (int elNumber = 0; elNumber < elementsNumber; elNumber++) {
		for (int i = 0; i < 4; i++) {
			grid.elements[elNumber].Ct0_dTplusP[i] = 0.0;
			for (int j = 0; j < 4; j++) {
				if (j == 0) {
					grid.elements[elNumber].Ct0_dTplusP[i] += (grid.elements[elNumber].C[i][j] / deltaTau) * t[grid.elements[elNumber].id[2] - 1];
				}
				if (j == 1) {
					grid.elements[elNumber].Ct0_dTplusP[i] += (grid.elements[elNumber].C[i][j] / deltaTau) * t[grid.elements[elNumber].id[3] - 1];
				}
				if (j == 2) {
					grid.elements[elNumber].Ct0_dTplusP[i] += (grid.elements[elNumber].C[i][j] / deltaTau) * t[grid.elements[elNumber].id[0] - 1];
				}
				

				if (j == 3) {
					grid.elements[elNumber].Ct0_dTplusP[i] += (grid.elements[elNumber].C[i][j] / deltaTau) * t[grid.elements[elNumber].id[1] - 1];
					grid.elements[elNumber].Ct0_dTplusP[i] += grid.elements[elNumber].P[i];
				}
			}
		}
	}
}

void SoE::printMatrixCt0_dTplusP(const Grid& grid, int elementsNumber)
{
	for (int elNumber = 0; elNumber < elementsNumber; elNumber++) {
		cout << "\nMatrix " << elNumber + 1 << " (CG/dT)*{t0} + {P}:" << endl;
		for (int i = 0; i < 4; i++) {
			cout << grid.elements[elNumber].Ct0_dTplusP[i] << "   ";
		}
		cout << endl;
	}
}

void SoE::zeroAggregatedMatrixP()
{
	for (int i = 0; i < n; i++) {
		PG[i] = 0.0;
	}
}

void SoE::aggregateMatrixCt0_dTplusP(const Grid& grid, int elementsNumber)
{
	for (int elementNumber = 0; elementNumber < elementsNumber; elementNumber++) {
		PG[grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].Ct0_dTplusP[0];
		PG[grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].Ct0_dTplusP[1];
		PG[grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].Ct0_dTplusP[2];
		PG[grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].Ct0_dTplusP[3];
	}
}

void SoE::initialStartTemperature(int t0)
{
	for (int i = 0; i < n; i++) {
		t[i] = t0;
	}
}

void SoE::solveSoE()
{
	// System of equations using the Gaussian elimination method
	
	// Forward Elimination (Eliminacja do przodu)
	for (int i = 0; i < n - 1; i++) {
		if (PG[i] == 0) {
			cout << "Error: Value '0' is on the diagonal of the matrix!" << endl;
			return;
		}

		for (int j = i + 1; j < n; j++) {
			double factor = HG[j][i] / HG[i][i]; //iloraz elementu macierzy w aktualnym wierszu i kolumnie przez element diagonalny
			for (int k = i; k < n; k++) {
				HG[j][k] -= factor * HG[i][k]; //eliminacja Gaussa poprzez odejmowanie wielokrotnoœci wiersza diagonalnego od aktualnego wiersza
			}
			PG[j] -= factor * PG[i]; //odejmowanie odpowiedniej wielokrotnoœci elementu P (wektora prawych stron) dla aktualnego wiersza
		}
	}

	// Back Substitution (Substytucja wsteczna)
	for (int i = n - 1; i >= 0; i--) {
		double sum = 0.0;
		for (int j = i + 1; j < n; j++) {
			sum += HG[i][j] * t[j]; //oblicza sumê iloczynów wspó³czynników elementów powy¿ej diagonalnej i odpowiadaj¹cych wartoœci wektora rozwi¹zania t
		}
		t[i] = (PG[i] - sum) / HG[i][i]; //rozwi¹zuje dla zmiennej w aktualnym wierszu przy u¿yciu wzoru substytucji wstecznej
	}
}

void SoE::printSoE()
{
	cout << "\nTemperature {t} n x 1:" << endl;
	for (int i = 0; i < n; i++) {
		cout << t[i] << "  ";
	}
	cout << endl;
}

void SoE::displayMinMaxTemperature(int time)
{
	double maxTemp = t[0];
	double minTemp = t[0];

	for (int i = 0; i < n; i++) {
		if (t[i] > maxTemp)
			maxTemp = t[i];
		if (t[i] < minTemp)
			minTemp = t[i];
	}
	cout << "\nTime = " << time << "s, minTemp = " << minTemp << ", maxTemp = " << maxTemp << endl;
}