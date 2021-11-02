#include "time.h"
#include <iostream>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include "Sim.cpp"

int main() {
	Sim *mySim = new Sim();
	mySim->MDSim();
	
	//delete [] mySim;
	puts("SIMULATION COMPLETE");
	
	
	return 0;
}
