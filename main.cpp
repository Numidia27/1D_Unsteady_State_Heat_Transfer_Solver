//#include "implicitSchemes.h"
#include "laasonen.h"
#include <iostream>
using namespace std;

int main()
{
	// Problem Variables (diffusivity, deltaX, deltaT, length, time, surfacetemp, initialtemp);
	Laasonen heatEquation(0.000217, 0.001, 0.01, 1.0, 0.18, 40.0, 0.0);

	// generate the .csv results file for laasonen scheme
	heatEquation.printLaasonen();

	cin.clear();
	return 0;
}
