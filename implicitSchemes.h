#ifndef IMPLICITSCHEMES_H // include guard
#define IMPLICITSCHEMES_H
#include <iostream>
#include <vector>
#include "parameters.h"
// #include <fstream>
using namespace std;

class ImplicitSchemes {
	protected:
		// The heat equation parameters
		double diffusivity, deltaX, deltaT, length, surfaceTemp, initialTemp, totalTime;
		int timePoints, spacePoints;

	public:
		// Default contructor
		ImplicitSchemes();
		
		// Alternate Constructor with arguments
		ImplicitSchemes(double diffusivity, double deltaX, double deltaT, double length, 
						double totalTime, double surfaceTemp, double initialTemp);

		// Getters and Setters of the parameters
		void setDiffusivity(double diffusivity);
		void setDeltaX(double deltaX);
		void setDeltaT(double deltaT);
		void setSpacePoints(int spacePoints);
		void setTimePoints(int timePoints);
		void setSurfaceTemp(double surfaceTemp);
		void setInitialTemp(double initialTemp);
		void setLength(double length);
		void setTotalTime(double totalTime);

		// Thomas Algorithm for solving the tridiagonal system A * x = d
		vector<double> thomasAlgorithm(vector<double> topDiagonal, vector<double> midDiagonal, 
										vector<double> botDiagonal, vector<double> d);
};
#endif
