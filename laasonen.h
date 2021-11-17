#ifndef LAASONEN_H // include guard
#define LAASONEN_H
#include <cmath>
#include "implicitSchemes.h"

class Laasonen : public ImplicitSchemes {
	protected:
		// Vectors of the tridiagonal matrix upper, middle and lower diagonals respectively
		vector<double>  upperDiagonal, mainDiagonal,lowerDiagonal;

		// Right hand side vector, that takes the boundary condition
		vector<double> temperature;
		
	public:
		// Default constructor
        Laasonen();

		// Alternate constructor
        Laasonen(double diffusivity, double deltaX, double deltaT, double length, 
				double totalTime, double surfaceTemp, double initialTemp);

		// Set Laasonen scheme parameters
		Laasonen setLaasonenParameters(Laasonen laasonen); 
		
		// Implement Laasonen scheme
		vector<double> laasonen();

		// Print the result in a .csv format
		void printLaasonen();
};
#endif
