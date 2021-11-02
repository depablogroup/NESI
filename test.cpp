#include "time.h"
#include <iostream>
#include <math.h>
#include <chrono>
#include <cstdlib>
#include <stdio.h>
#include "Sim.cpp"

using namespace NESI;
using namespace std;

int main() {
	Sim *mySim = new Sim();

    auto start = std::chrono::high_resolution_clock::now(); 

    for (int i = 0; i < 1000; i++ ) {
        mySim->calc_nonbond_forces_omp(true,true);
    }
    
    auto stop = std::chrono::high_resolution_clock::now();
    std:chrono::duration<double> duration = stop - start;
    printf("Average execution time: %lf ms\n",duration.count()); 
    
    return 0;
}
