#include "math_vector.h"

using namespace NESI;

class Box {
	public:
            
        // Box dimensions
        double x;
        double y;
        double z;
        int hasWallZ;
        double slabWidth;
        int numSlabs;
        int numInsertionsPerStep;

        // some system data
        double tempSet; // desired temperature
        double tempSetTop; // desired temperature at the top wall
        double velVar; // variance of velocity distribution
        double velVarTop; // variance of velocity distribution at the top wall
        double targetK;  // target kinetic energy
        double rescaleTime; // rescaling time for lower thermostat (CSVR thermostat)
        double rescaleTimeTop;  // rescale time for top thermostat (Andersen thermostat)
        double collisionProb; // collision probability for andersen thermostat (upper)
        double thermoHeight; // thermostat/chem potential bath height at the bottom
        double thermoHeightTop; // thermostat/chem potential bath height at bottom
        int zeroLiqMomentum;
        int doThermostat;
        double pinningForce;
        double pinHeight;
        double wallEps;
        double wallCut;
        double temp; // actual (instantaneous) temperature
        double dens;  //Inverse density of the box
        double E_tot; //Total energy
        int seed; //Random number seed
        int numBeads; //Number of atoms in the simulation
        int numTypes;     //Number of different types of interacting particles in the simulation  
        double deltaMuTop;  // target chemical potential difference at top of system
        double deltaMuBottom; // target chemical potential difference at bottom of system
        int swapEvery; // how often to try to swap particles in the baths
        int upperAccept;
        int lowerAccept;

        // linked-list data
        double rcut; //Cutoff distance for energy calculations and grid/cell width
        int cellMult;      // number of cells in each cutoff distance
        int cell_lim_x;        // how many cells in x direction are required for force calculations
        int cell_lim_y;        // how many cells in y direction are required for force calculations
        int cell_lim_z;        // how many cells in z direction are required for force calculations
        double cell_width_x; //These things specify how far we have to search with each grid cell when computing the neighbor list
        double cell_width_y; // must be at least rcut
        double cell_width_z;
        int cell_num_x;    //These three are determined when the simulation is initialized
        int cell_num_y;
        int cell_num_z;
        int numCell;       //Number of total grid sites in the simulation
        int *head; // array of list heads, one for each cell. The value of head[cellNumber] is the first atom in the cell's list below
        int *bead_list; // linked list of atoms in cells - atom_list[i] is the index of the next atom in the list
        lvector *cellVec; // x,y,z indices of the cell
        int **cell_neighbors;  // array of arrays that gives neighbors for each cell
        bool **do_force_calc;  // for each cell in the list of neighbors, only half need to be included in force calculation (Newton's law)
        int *cell_num_neighbors; // array giving number of neighbors for each cell (may differ depending on PBCs)
        int maxNeighbors;       // max neighbors possible for a cell

            
        // run control, set in input file
        int numSteps; //Number of MD steps
        int numEqSteps; //Number of equilibration steps
        double dt; // time step in LJ units
        int step; //Current step
        int calcInterval; //Interval after which data is calculated
        int writeInterval; // interval after which data is written/averaged
        int writeCoords;  // interval after which coordinates are written
        int checkpoint; //0 if starting from scratch. 1 if loading a configuration
        int checkInterval; //MD step interval after which checkpoints are made
        int giveNewVel; // given new velocities
        int continuation;

        inline bool doNeighborCalc(lvector &a, lvector &b) {
            
            if (a[2] < b[2]) return false;
            if (a[2] == b[2] && a[1] < b[1]) return false;
            if (a[2] == b[2] && a[1] == b[1] && a[0] < b[0]) return false;
            return true;
        }

};
