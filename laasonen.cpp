#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "laasonen.h"

using namespace std;

//CONSTRUCTORS
/*
* Default constructor - that inherits from the ImplicitSchemes default constractor, @see ImplicitSchemes
*/
Laasonen::Laasonen() {
	lowerDiagonal.push_back(0);
	mainDiagonal.push_back(1);
	upperDiagonal.push_back(0);
}

/*
* Alternate constructor - that inherits from the ImplicitSchemes alternate constructor, @see ImplicitSchemes
*/ 
Laasonen::Laasonen(double diffusivity, double deltaX, double deltaT, double length, 
					double totalTime, double surfaceTemp, double initialTemp) 
					: ImplicitSchemes(diffusivity, deltaX, deltaT, length, totalTime, surfaceTemp, initialTemp) {}

/*
* Initialize and set all the parameters of Laasonen scheme
*/
Laasonen Laasonen::setLaasonenParameters(Laasonen laasonen) {
	laasonen.setDiffusivity(diffusivity);
	laasonen.setDeltaX(deltaX);
	laasonen.setDeltaT(deltaT);
	laasonen.setSpacePoints(int(length / deltaX)); 
	laasonen.setTimePoints(int(totalTime / deltaT)); 
	laasonen.setSurfaceTemp(surfaceTemp);
	laasonen.setInitialTemp(initialTemp);
	return laasonen;
}

/*
* Implement Laasonen Scheme to solve the 1D Heat conduction equation
*/
vector<double> Laasonen::laasonen() {
	/*
	* brief solve 1D heat equation using Laasonen leads to a tridiagonal system 
	* that can be efficiently solve by Thomas Algorithm, @see the implementation of Thomas Algorithm in ImplicitSchemes file
	*/

	// @param a is the constant of the 1D heat partial differential equation
	double a = diffusivity * (deltaT / pow(deltaX, 2));

	/* 
	* Initialize the first element of the right hand side vector d, of the system A*x = d
	* to the left side surface temperature of the wall, which is part of the boundary condition
	*/
	temperature.push_back(surfaceTemp);

	/*
	* Setting the value of the diagonals of the Laasonen matrix system 
	* as explained in the report, the upper and lower diagonals set to (-a), whereas the main diagonal set to (1 + (2 * a))
	*/
	for (int i = 1; i < spacePoints; i++) {
		lowerDiagonal.push_back(-a);
		mainDiagonal.push_back(1 + (2 * a));
		upperDiagonal.push_back(-a);
		temperature.push_back(initialTemp);
	}

	// lowerDiagonal[0] doesn't exist
	// upper[0] = 0
	upperDiagonal.push_back(0); 

	// mainDiagonal[0] = 1
	mainDiagonal.push_back(1);	 

	/* 
	* Initialize the last element of the right hand side vector d, of the system A*x = d
	* to the right side surface temperature of the wall, which is part of the boundary condition
	*/
	temperature.push_back(surfaceTemp);

	/* 
	* Solve Laasonen scheme tridiagonal system using Thomas Algorithm
	* @return temperature
	*/
	for(int t = 0; t < timePoints; t++) {
		temperature = thomasAlgorithm(lowerDiagonal, mainDiagonal, upperDiagonal, temperature);
	}
	return temperature;
}

// PRINT RESULTS
/*
* Print the results of Laasonen scheme in .csv file
*/
void Laasonen::printLaasonen() {
	Laasonen laasonenPrintResult;
	laasonenPrintResult = Laasonen::setLaasonenParameters(laasonenPrintResult);
	vector<double> vectResult = laasonenPrintResult.laasonen();
	ofstream resultsFile("laasonen.csv");
	if (resultsFile.is_open()) {
	
		for (unsigned i = 0; i < vectResult.size(); i++) {
			resultsFile << (i * deltaX) << "	" << vectResult[i] << endl;
		}
		resultsFile.close();
	}
}
