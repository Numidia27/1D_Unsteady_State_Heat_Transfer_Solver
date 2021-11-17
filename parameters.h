#ifndef PARAMETERS_H // include guard
#define PARAMETERS_H

class Parameters
{
	double diffusivity = 0.1; // 0.1 ft^2/hr
	double deltaX = 0.05; // 0.05
	double deltaT = 0.01; //0.01
	//vector<double> vecDeltaT
	double length = 1; // 1 ft
	int spacePoints = int((length / deltaX) + 0.5); // 40
	//vector<int> timePoints; // 0, 0.5
	int timePoints = 50;
	double surfaceTemp = 300; // 300
	double initialTemp = 100; // 100
};
#endif
