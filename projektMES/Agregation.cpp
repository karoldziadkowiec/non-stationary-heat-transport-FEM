#include <cmath>

#include "Grid.h"
#include "Agregation.h"

SoE::SoE(int n)
{
	this->n = n;

	HG = new double* [n];
	for (int i = 0; i < n; i++) {
		HG[i] = new double[n] {};
	}
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			HG[i][j] = 0.0;
		}
	}

	P = new double[n];
	t = new double[n];
	for (int i = 0; i < n; i++) {
		P[i] = 0.0;
		t[i] = 0.0;
	}
}

SoE::~SoE()
{
	for (int i = 0; i < n; i++) {
		delete[] HG[i];
	}
	delete[] HG;
	delete[] P;
}

void SoE::aggregateMatrixH(const Grid& grid, int elementNumber)
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

void SoE::printAggregatedMatrixH()
{
	cout << "\nAgrregated Matrix [HG] n x n:" << endl;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			cout << HG[i][j] << " ";
		}
		cout << endl;
	}
}

void SoE::aggregateVectorP(const Grid& grid, int elementNumber)
{
	P[grid.elements[elementNumber].id[2] - 1] += grid.elements[elementNumber].P[0];
	P[grid.elements[elementNumber].id[3] - 1] += grid.elements[elementNumber].P[1];
	P[grid.elements[elementNumber].id[0] - 1] += grid.elements[elementNumber].P[2];
	P[grid.elements[elementNumber].id[1] - 1] += grid.elements[elementNumber].P[3];
}

void SoE::printAggregatedVectorP()
{
	cout << "\nAgrregated Vector {P} n x 1:" << endl;
	for (int i = 0; i < n; i++) {
		cout << P[i] << "  ";
	}
	cout << endl;
}

void SoE::solveSoE()
{
	// System of equations using the Gaussian elimination method
	
	// Forward Elimination (Eliminacja do przodu)
	for (int i = 0; i < n - 1; i++) {
		if (P[i] == 0) {
			cout << "Error: Value '0' is on the diagonal of the matrix!" << endl;
			return;
		}

		for (int j = i + 1; j < n; j++) {
			double factor = HG[j][i] / HG[i][i]; //iloraz elementu macierzy w aktualnym wierszu i kolumnie przez element diagonalny.
			for (int k = i; k < n; k++) {
				HG[j][k] -= factor * HG[i][k]; //eliminacja Gaussa poprzez odejmowanie wielokrotnoœci wiersza diagonalnego od aktualnego wiersza.
			}
			P[j] -= factor * P[i]; //odejmowanie odpowiedniej wielokrotnoœci elementu P (wektora prawych stron) dla aktualnego wiersza.
		}
	}

	// Back Substitution (Substytucja wsteczna)
	for (int i = n - 1; i >= 0; i--) {
		double sum = 0.0;
		for (int j = i + 1; j < n; j++) {
			sum += HG[i][j] * t[j]; //oblicza sumê iloczynów wspó³czynników elementów powy¿ej diagonalnej i odpowiadaj¹cych wartoœci wektora rozwi¹zania t.
		}
		t[i] = (P[i] - sum) / HG[i][i]; //rozwi¹zuje dla zmiennej w aktualnym wierszu przy u¿yciu wzoru substytucji wstecznej.
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