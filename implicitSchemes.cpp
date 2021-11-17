#include "implicitSchemes.h"
#include <cmath>

using namespace std;

// CONSTRUCTORS
/*
* Default constructor
*/
ImplicitSchemes::ImplicitSchemes() {}

/*
* Alternate constructor - creates a vector and assign values to the parameters
*/
ImplicitSchemes::ImplicitSchemes(double diffusivity, double deltaX, double deltaT, 
double length, double totalTime, double surfaceTemp, double initialTemp){
	this->diffusivity = diffusivity;    
	this->deltaX = deltaX;             
	this->deltaT = deltaT;             
	this->length = length;              
	this->totalTime = totalTime;    
	this->surfaceTemp = surfaceTemp;         
	this->initialTemp = initialTemp;	 
}

/*
* Accessor methods - set the parameters of the heat equation
*/
void ImplicitSchemes::setDiffusivity(double diffusivity) {
	this->diffusivity = diffusivity;
}

void ImplicitSchemes::setDeltaX(double deltaX) {
	this->deltaX = deltaX;
}

void ImplicitSchemes::setDeltaT(double deltaT) {
	this->deltaT = deltaT;
}

void ImplicitSchemes::setSpacePoints(int spacePoints) {
	this->spacePoints = spacePoints;
}

void ImplicitSchemes::setTimePoints(int timePoints) {
	this->timePoints = timePoints;
}

void ImplicitSchemes::setSurfaceTemp(double surfaceTemp) {
	this->surfaceTemp = surfaceTemp;
}

void ImplicitSchemes::setInitialTemp(double initialTemp) {
	this->initialTemp = initialTemp;
}

void ImplicitSchemes::setLength(double length) {
	this->length = length;
}

void ImplicitSchemes::setTotalTime(double totalTime) {
	this->totalTime = totalTime;
}

// SOLVE TRIDIAGONAL SYSTEM A * x = d
/*
* Implement Thomas Algorithm to solve the tri-diagonal system A * x = d of the implicit schemes
*/
vector<double> ImplicitSchemes::thomasAlgorithm(vector<double> topDiagonal, vector<double> midDiagonal, 
												vector<double> botDiagonal, vector<double> d) {												
	/* 
	* @param n is the size in the tri-diagonal system A * x = d
	* n is the size of the vector midDiagonal and the right hand side vector d
	* n-1 is the size of the vector topDiagonal and botDiagonal
	*/
	int n = d.size(); 
	
	// x is the vector solution of the unknowns to calculate in the tri-diagonal system A * x = d
	vector<double> x(n,0);

	/*
	* brief Thomas Algorithm based on LU decomposition, where the system A*x = d can be re-written as LU = r
	* where L, U are the lower and upper triangular matrices respectively
	* then the system can be solved, first by solving Lp = r for p and then Ux = p for x, which is the solution
	* Forwards steps:
	* Solve the Lp = r, by setting the main diagonal to 1 and lower and upper diagonals to 0
	*/
	for (int i = 0; i < n; i++) {
		if (i == 0){
			topDiagonal[0] /= midDiagonal[0];
			d[0] /= midDiagonal[0];
			midDiagonal[0] /= midDiagonal[0]; 
		} else {
			midDiagonal[i] -= (topDiagonal[i - 1] * botDiagonal[i]); 
			d[i] -= (d[i - 1] * botDiagonal[i]);
			botDiagonal[i] = 0; 
			topDiagonal[i] /= midDiagonal[i];
			d[i] /= midDiagonal[i];
		}	
	}
	/* 
	* Backwards steps:
	* Solve the Ux = p for x
	* @return x, the solution of A * x = d
	*/
	for (int i = n-1; i >= 0; i--) {
		if (i == d[n]) {
			x[n] = d[n];
		}
		else {
			x[i] = d[i] - (topDiagonal[i] * x[i + 1]);
		}
	}
	return x;
}
