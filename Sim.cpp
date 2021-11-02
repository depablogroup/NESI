#include "time.h"
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <cstdlib>
#include <stdio.h>
#include <chrono>
#include <algorithm>
#include "Sim.hpp"
#include "mersenne.h"
#include <omp.h>
#include <climits>
#define PI 3.14159265358979323846
#define xi1 0.099194
#define xi2 -0.16346
#define rs 1.791044776 
#define EMPTY -1
#define xi1Wall -0.58202813
#define xi2Wall -1.38057246

using namespace NESI;
using namespace std;

/****************************** SIMULATION RUN CONTROL ********************************/


// system preparation and initialization
Sim::Sim() {

    read_input(); // get input paramters
    if (box.checkpoint) read_configuration();
    else init_configuration(); // get starting configuration
    wrap_atoms();
    generate_matrices(); // set up interactions
    initialize_grid(); // build cell lists
    initialize_data(); // get ready to store data
    if (box.giveNewVel) set_velocities(); // set velocities
    if (!box.continuation) openFiles();  // prepare to write output
}

// perform simulation
void Sim::MDSim() {
   
    // start with a force calculation
    calcVirials = calcEnergy = true;
    calc_forces(calcVirials,calcEnergy);
    box.upperAccept = box.lowerAccept = 0;

    // start the timer
    auto start = std::chrono::high_resolution_clock::now();

    
    //equilibration - slab data is not calculated here, only coordinates are output
    for (box.step=0;box.step<box.numEqSteps;box.step++) {

        // do a half-step of velocity update
        update_vel_half(vs,fs,box.dt);

        // do a full step of position update
        update_pos(xs,vs,box.dt);

        // rebuild the lists and calculate forces
	build_bead_lists();
        calcVirials = false;
        if (box.swapEvery > 0) calcEnergy = ((box.step+1)%box.swapEvery == 0);
        calc_forces(calcVirials,calcEnergy);
	
        // final half-step velocity update
        update_vel_half(vs,fs,box.dt);
     
        // thermostat the portion near the lower wall
        if (box.hasWallZ) scale_wall_velocities();
        else scale_velocities();
        
        // write coordinates
        if ((box.step+1)%box.writeCoords == 0) {
            printf("Turn: %d\n",box.step+1);
            printXYZ();
        }
        
        // particle swaps to maintain chemical potential differences
        if ((box.swapEvery > 0) && ((box.step+1)%box.swapEvery == 0)) swap_particles();
        
        // zero out the liquid momentum if necessary
        if (box.zeroLiqMomentum > 0 && (box.step+1)%box.zeroLiqMomentum == 0) zero_wall_momentum();

    }

    // production
    for (box.step=box.numEqSteps;box.step<box.numSteps+box.numEqSteps;box.step++) {
        
        // do a half-step of velocity update
        update_vel_half(vs,fs,box.dt);

        // do a full step of position update
        update_pos(xs,vs,box.dt);
        
        // rebuild the lists and calculate forces
        calcVirials = ((box.step+1)%box.calcInterval == 0);
        calcEnergy = (((box.step+1)%box.calcInterval == 0) || ((box.step+1)%box.swapEvery == 0));
        build_bead_lists();
        calc_forces(calcVirials,calcEnergy);

        // final half-step velocity update
        update_vel_half(vs,fs,box.dt);

        // apply thermostat at walls
        if (box.hasWallZ) scale_wall_velocities();
        else scale_velocities();

       
        // calculate slab data
        if ((box.step+1)%box.calcInterval == 0) calc_slab_data_omp();

        // write slab data
        if ((box.step+1)%box.writeInterval == 0) {
            write_slab_data();
        }

        // write coordinates
        if ((box.step+1)%box.writeCoords == 0) {
            printf("Turn: %d\n",box.step+1);
            printXYZ();
        }

        // write checkpoint file and widom histograms
        if ((box.step+1)%box.checkInterval == 0) {
            print_configuration();
        }

        // particle swaps to maintain chemical potential differences
	// this messes with the per-particle energies, so must be done AFTER the slab calculations
        if ((box.swapEvery > 0) && ((box.step+1)%box.swapEvery == 0)) swap_particles();

        // zero out the liquid momentum if necessary
        if (box.zeroLiqMomentum > 0 && (box.step+1)%box.zeroLiqMomentum == 0) zero_wall_momentum();


    }
    
    cout << "Swaps accepted at upper bath: " << box.upperAccept << endl;
    cout << "Swaps accepted at lower bath: " << box.lowerAccept << endl;


    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    double tsps = (box.numSteps + box.numEqSteps) / (duration.count());
    printf("runtime %f\n%e timesteps per second\n",duration.count(), tsps);

    //write_slab_data();
    dumpData();
    //printXYZ();
}

/******************************************************************************/


/********************** SYSTEM PREPARATION AND INPUT *********************************/

// prepare the files we're going to write to
void Sim::openFiles() {
    FILE *dump;
    dump = fopen("config.xyz","w");
    fclose(dump);
    
    dump = fopen("data.txt","w");
    fclose(dump);

    int i;
    for (i=0;i<box.numTypes;i++) {
        char filename[32];
        sprintf(filename,"slab_dens_%d.txt",i);
        dump=fopen(filename,"w");
        fclose(dump);
    }
    for (i=0;i<box.numTypes;i++) {
        char filename[32];
        sprintf(filename,"slab_mtms_%d.txt",i);
        dump=fopen(filename,"w");
        fclose(dump);
    }
    dump=fopen("slab_kins.txt","w");
    fclose(dump);
    dump=fopen("slab_pots.txt","w");
    fclose(dump);
    dump=fopen("slab_virs.txt","w");
    fclose(dump);
    dump=fopen("slab_surf.txt","w");
    fclose(dump);
    dump=fopen("slab_heat_fluxes.txt","w");
    fclose(dump);

    for (i=0;i<box.numTypes;i++) {
        char filename[32];
        sprintf(filename,"slab_widom_hist_%d.txt",i);
        dump=fopen(filename,"w");
        fclose(dump);
    }

}

void Sim::read_input() {

  FILE *input;
  char tt[2001];
  int ret;
  input=fopen("input", "r");

  if (input!=NULL)
    { 	/******		Reading Input File		********/
      ret=fscanf(input,"%d", &box.seed);		     fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.checkpoint);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.continuation);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.giveNewVel);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numSteps);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numEqSteps);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.calcInterval);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.writeInterval);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.writeCoords);		 fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.checkInterval);    fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.x);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.y);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.z);               fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numBeads);               fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numTypes);               fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.cellMult);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.tempSet);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.tempSetTop);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.thermoHeight);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.thermoHeightTop);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.dt);               fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.hasWallZ);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.slabWidth);               fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.numInsertionsPerStep);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &minEnergy);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &maxEnergy);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &binWidth);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.rescaleTime);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.rescaleTimeTop);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.pinningForce);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.pinHeight);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.wallEps);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.wallCut);               fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.zeroLiqMomentum);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.deltaMuTop);               fgets(tt,2000,input);
      ret=fscanf(input,"%lf", &box.deltaMuBottom);               fgets(tt,2000,input);
      ret=fscanf(input,"%d", &box.swapEvery);               fgets(tt,2000,input);
      fgets(tt,2000,input);fgets(tt,2000,input);
      box.velVar = sqrtl(box.tempSet);
      box.velVarTop = sqrtl(box.tempSetTop);
      box.numSlabs = (int) floor(box.z/box.slabWidth);
      box.slabWidth = box.z/box.numSlabs;
      box.targetK = 1.5*box.numBeads*box.tempSet;
      box.collisionProb = box.dt/box.rescaleTimeTop;
      // set histogram parameters
      quitEnergy = maxEnergy-minEnergy;
      numBins = (int) ceil(quitEnergy/binWidth); // number of actual bins is one larger to handle all high-energy values
      binWidth = quitEnergy/numBins;

      fclose(input);
    }

  else
    { 
      fprintf(stdout,"\n Input File Not Found\n");
      exit(EXIT_FAILURE);
    }
   
    // get openMP parameters
    #if defined(_OPENMP)
    #pragma omp parallel 
    {
        if (omp_get_thread_num() == 0) {
            nThreads = omp_get_num_threads();
            cout << "Number of threads: " << nThreads << endl;
        }
    }
    #else
    nThreads = 1;
    #endif


}

// read configuration from initial state (provided as input)
void Sim::read_configuration() {
    int i,start,end,anchor,rigid,issite,hassites;
    char tt[2001];
    double total_volume = 0.0;

    // xyz-input file with velocities
    FILE *atom_file = fopen("init.xyz","r"); //Location of all atom files 
    //fgets(tt,2000,atom_file);             //Section title
    fscanf(atom_file, "%d",&box.numBeads);
    fscanf(atom_file, "%lf%lf%lf", &box.x, &box.y, &box.z); // comment line is xyz dimensions
    xs = new NESI::vector[box.numBeads];
    vs = new NESI::vector[box.numBeads];
    fs = new NESI::vector[box.numBeads*nThreads];
    types = new int[box.numBeads];
    virials = new shape[box.numBeads*nThreads];
    potentials = new double[box.numBeads*nThreads];
    
    // read init config file
    for (i=0;i<box.numBeads;i++) {
        fscanf(atom_file,"%d%lf%lf%lf%lf%lf%lf", &types[i], &xs[i][0], &xs[i][1], &xs[i][2], &vs[i][0], &vs[i][1], &vs[i][2]);
    }
    fclose(atom_file);
    
    wrap_atoms();
    double volume = box.x*box.y*box.z;
    box.dens = box.numBeads/volume;


}

// initialize configuration by placing atoms on a lattice
void Sim::init_configuration() {
    int i,j,k,dim,idx,factor,dimZ;
    double dx,dy,dz;
    double spacing = 1.2;
    char tt[2001];


    // initialize the data structures
    types = new int[box.numBeads];
    xs = new NESI::vector[box.numBeads];
    vs = new NESI::vector[box.numBeads];
    fs = new NESI::vector[box.numBeads*nThreads];
    virials = new shape[box.numBeads*nThreads];
    potentials = new double[box.numBeads*nThreads];

    // find lattice spacing needed
    factor = (int) box.z/box.x;
    dim = (int) ceil(cbrt(box.numBeads/factor));  // now many beads in the x and y directions
    dx = box.x/dim;
    dy = box.y/dim;
    dz = box.z/(dim*factor);

    // loop through atoms, placing them on the lattice
    for (i=0;i<dim;i++) {
        for (j=0;j<dim;j++) {
            for (k=0;k<factor*dim;k++) {
                
                // exit if we have enough atoms
                idx = i*dim*dim*factor + j*dim*factor + k;
                if (idx >= box.numBeads) continue;

                // place the bead at the center of the lattice
                types[idx] = 0; 
                xs[idx][0] = ((double)i+0.5)*dx;
                xs[idx][1] = ((double)j+0.5)*dy; 
                xs[idx][2] = ((double)k+0.5)*dz; 

            }
        }
    }
  
    wrap_atoms();
    double volume = box.x*box.y*box.z;
    box.dens = box.numBeads/volume;
}

// draw velocities at random from boltzmann dist.
void Sim::set_velocities() {
    
    int i,j;
    double var = sqrtl(box.tempSet);
    double kinetic = 0.0;
    double temp,scale;
    NESI::vector p,pTot;

    // generate random velocities and keep track of total momentum - assumes particles of identical mass
    for (i=0;i<box.numBeads;i++) {
        for (j=0;j<3;j++) {
            vs[i][j] = calc_random_gaussian(0.0,var*sqrtInvMasses[types[i]]);
        }
        vec_scalar_mult(p,vs[i],masses[types[i]]);
        vec_add(pTot,p);
    }
    vec_scalar_mult(pTot,1.0/((double) box.numBeads));
    
    // subtract off the total momentum and keep track of kinetic energy
    for (i=0;i<box.numBeads;i++){
        vec_scalar_mult(p,pTot,invMasses[types[i]]);
        vec_subtr(vs[i],p);
        kinetic += 0.5*masses[types[i]]*vec_dot(vs[i],vs[i]);
    }

    // rescale velocities to get correct temp
    scale = sqrtl(box.targetK/kinetic);
    for (i=0;i<box.numBeads;i++){
        vec_scalar_mult(vs[i],scale);
    }

}

// read in the interaction data and create parameter matrices
void Sim::generate_matrices() {
    
    int ret,i,j,type1,type2;
    double epsVal,sigVal,rcutVal,massVal;
    const char* chr;
    
    // set up the data structures
    eps = new double *[box.numTypes];
    sigs = new double *[box.numTypes];
    rcuts = new double *[box.numTypes];
    masses = new double[box.numTypes];
    invMasses = new double[box.numTypes];
    sqrtInvMasses = new double[box.numTypes];

    for (i=0;i<box.numTypes;i++) {
        
        eps[i] = new double [box.numTypes];
        sigs[i] = new double [box.numTypes];
        rcuts[i] = new double [box.numTypes];

    }
 
    // set everything to zero initially
    for (i=0;i<box.numTypes;i++) {
        for (j=0;j<box.numTypes;j++) {
            eps[i][j] = 0.0;
            sigs[i][j] = 0.0;
            rcuts[i][j] = 0.0;
        }
        masses[i] = invMasses[i] = sqrtInvMasses[i] = 1.0;  // just in case, since we divide by this pretty often
    }


    // read in the mass data line-by-line
    std::ifstream typeData("types.dat");
    for (std::string line; getline(typeData,line);) {
        
        // convert string to char*
        chr = line.c_str();

        // read in the LJ types the parameters apply to
        ret=sscanf(chr,"%d%lf", &type1,&massVal);		     
        masses[type1] = massVal; 
        invMasses[type1] = 1.0/massVal;
        sqrtInvMasses[type1] = sqrtl(invMasses[type1]);

    }

    // read in the interaction data line-by-line
    std::ifstream data("interactions.dat");
    for (std::string line; getline(data,line);) {
        
        // convert string to char*
        chr = line.c_str();

        // read in the LJ types the parameters apply to
        ret=sscanf(chr,"%d%d%lf%lf%lf", &type1,&type2,&epsVal,&sigVal,&rcutVal);		     
        
        eps[type1][type2] = eps[type2][type1] = epsVal;
        sigs[type1][type2] = sigs[type2][type1] = sigVal;
        rcuts[type1][type2] = rcuts[type2][type1] = rcutVal;
        
    }
    
    data.close();

    // use appropriate combination rules for any missing elements
    for (i=0;i<box.numTypes;i++) {
        for(j=i+1;j<box.numTypes;j++) {
       
            if (eps[i][j] == 0.0) eps[i][j] = eps[j][i] = sqrtl(eps[i][i]*eps[j][j]);
            if (sigs[i][j] == 0.0 && sigs[i][i] != 0.0 && sigs[j][j] != 0.0) sigs[i][j] = sigs[j][i] = 0.5*(sigs[i][i]+sigs[j][j]);
            if (rcuts[i][j] == 0.0 && rcuts[i][i] != 0.0 && rcuts[j][j] != 0.0) rcuts[i][j] = rcuts[j][i] = max(rcuts[i][i],rcuts[j][j]);

        }
    }

    // find the minimum value of epsilon (used in exit criteria for widom insertion)
    minEps = 9999999999.9;
    for (i=0;i<box.numTypes;i++) {
        for (j=0;j<box.numTypes;j++) {
            if (eps[i][j] < minEps) minEps = eps[i][j];
        }
    }
    

}

// set up random num. generators and slab data structures
void Sim::initialize_data() {
    
    if (box.seed < 0) {
        box.seed = time(NULL);
    }
    RanGen = new CRandomMersenne(box.seed);
    threadRanGen = (CRandomMersenne*) malloc(nThreads*sizeof(CRandomMersenne));
    int i,j;
    for (i=0;i<nThreads;i++) {
        threadRanGen[i] = *new CRandomMersenne(box.seed+i+1); // make sure each RNG has a different seed
    }


    // initialize slab data
    slab_counts = new int *[box.numTypes];
    for (i=0;i<box.numTypes;i++) {
        slab_counts[i] = new int [box.numSlabs*nThreads];
    }
    slab_mtms = new NESI::vector *[box.numTypes];
    for (i=0;i<box.numTypes;i++) {
        slab_mtms[i] = new NESI::vector [box.numSlabs*nThreads];
    }
    slab_kins = new shape[box.numSlabs*nThreads];
    slab_pots = new double[box.numSlabs*nThreads];
    slab_heat_fluxes = new NESI::vector[box.numSlabs*nThreads];
    slab_virs = new shape[box.numSlabs*nThreads];
    typeEnergies = new double[box.numTypes*nThreads]; 
    slab_widom_hist = new long int **[box.numSlabs];
    for (i=0;i<box.numTypes;i++) {
        slab_widom_hist[i] = new long int *[box.numSlabs];
        for (j=0;j<box.numSlabs;j++) {
            slab_widom_hist[i][j] = new long int [numBins+1]; // extra bin to handle super-high energy results that essentially just contribute zeros in the boltzmann factor
        } // note that each thread has its own chunk of slabs for insertions, so no write conflicts
    }

    // just to be safe, initialize all of these to zeros
    for (i=0;i<box.numTypes;i++) {
        memset(slab_counts[i], 0, box.numSlabs*nThreads*sizeof(*slab_counts[i]));
        memset(slab_mtms[i], 0, box.numSlabs*nThreads*sizeof(*slab_mtms[i]));
        for (j=0;j<box.numSlabs;j++) {
            memset(slab_widom_hist[i][j], 0, (numBins+1)*sizeof(*slab_widom_hist[i][j]));
        }
    }
    memset(slab_kins, 0, box.numSlabs*nThreads*sizeof(*slab_kins));
    memset(slab_pots, 0, box.numSlabs*nThreads*sizeof(*slab_pots));
    memset(slab_virs, 0, box.numSlabs*nThreads*sizeof(*slab_virs));
    memset(slab_heat_fluxes, 0, box.numSlabs*nThreads*sizeof(*slab_heat_fluxes));

}

/*********************************************************************************/


/********************* CELL LIST METHODS ******************************************/

// maps a position (vector) to a given cell
int Sim::mapToGrid(NESI::vector r) {
    double x,y,z,x1,y1,z1;
    int a,b,c,m;

    x = r[0];
    y = r[1];
    z = r[2];

    a = floor(x/box.cell_width_x);
    b = floor(y/box.cell_width_y);
    c = floor(z/box.cell_width_z);
    
    m = c*box.cell_num_x*box.cell_num_y + b*box.cell_num_x + a;

    return m;
}

// maps a position (x,y,z) to a given cell
int Sim::mapToGrid(double x, double y, double z) {
    int a,b,c,m;

    a = floor(x/box.cell_width_x);
    b = floor(y/box.cell_width_y);
    c = floor(z/box.cell_width_z);
    
    m = c*box.cell_num_x*box.cell_num_y + b*box.cell_num_x + a;

    return m;
}

// converts cell coordinates into position in cell data list
int Sim::calc_cell_placement(int a, int b, int c) {
    a = a%box.cell_num_x;
    if (a<0) a += box.cell_num_x;
    b = b%box.cell_num_y;
    if (b<0) b += box.cell_num_y;
    c = c%box.cell_num_z;
    if (c<0) c += box.cell_num_z;
    
    return c*box.cell_num_y*box.cell_num_x + b*box.cell_num_x + a;
}


// converts position in cell data list into cell coordinates
void Sim::calc_cell_index(NESI::lvector &index, int m) {
    int a,b,c;

    c = m/(box.cell_num_x*box.cell_num_y);
    m -= c*box.cell_num_x*box.cell_num_y;
    
    b = m/box.cell_num_x;
    a = m - b*box.cell_num_x;

    index[0] = a;
    index[1] = b;
    index[2] = c;
}

// create the cell list data structures
void Sim::initialize_grid() {

    int i,j,k,m1,m2;
    lvector idx1,idx2;

    // first get longest cut off distance
    box.rcut = 0.0;
    for (i=0;i<box.numTypes;i++) {
        for (j=i;j<box.numTypes;j++) {
            if (rcuts[i][j] >= box.rcut) box.rcut = rcuts[i][j];
        }
    }
   
    // figure out how many cells we have, how big they need to be, and how 
    // many neighbors are needed in force calculations
    box.cell_num_x = int(floor(box.x * box.cellMult / box.rcut));
    box.cell_num_y = int(floor(box.y * box.cellMult / box.rcut));
    box.cell_num_z = int(floor(box.z * box.cellMult / box.rcut));

    box.cell_width_x = box.x / double(box.cell_num_x);
    box.cell_width_y = box.y / double(box.cell_num_y);
    box.cell_width_z = box.z / double(box.cell_num_z);

    box.numCell = box.cell_num_x*box.cell_num_y*box.cell_num_z;
    box.cell_lim_x = int(ceil(box.rcut / box.cell_width_x));
    box.cell_lim_y = int(ceil(box.rcut / box.cell_width_y));
    box.cell_lim_z = int(ceil(box.rcut / box.cell_width_z));

    // construct the list of cell neighbors
    box.maxNeighbors = (2*box.cell_lim_x+1)*(2*box.cell_lim_y+1)*(2*box.cell_lim_z+1) - 1;
    box.cell_neighbors = new int *[box.numCell];
    box.cell_num_neighbors = new int[box.numCell];
    box.do_force_calc = new bool *[box.numCell];
    box.cellVec = new lvector[box.numCell];
    for (m1=0;m1<box.numCell;m1++) {
        
        // initialize list of neighbors
        box.cell_neighbors[m1] = new int[box.maxNeighbors];
        box.do_force_calc[m1] = new bool[box.maxNeighbors];
        box.cell_num_neighbors[m1] = 0;

        // get cell index 
        calc_cell_index(idx1,m1);
        box.cellVec[m1][0] = idx1[0];
        box.cellVec[m1][1] = idx1[1];
        box.cellVec[m1][2] = idx1[2];

        // loop through nearby cells to find neighbors
        for (i=idx1[0]-box.cell_lim_x;i<=idx1[0]+box.cell_lim_x;i++) {
            for (j=idx1[1]-box.cell_lim_y;j<=idx1[1]+box.cell_lim_y;j++) {
                for (k=idx1[2]-box.cell_lim_z;k<=idx1[2]+box.cell_lim_z;k++) {

                    // if we have walls in the z-direction, ignore this cell and move on
                    if (box.hasWallZ && (k < 0 || k >= box.cell_num_z)) continue;
                    if (k == idx1[2] && j == idx1[1] && i == idx1[0]) continue; // don't include ourselves in the list

                    // get list placement of the neighboring cell
                    m2 = calc_cell_placement(i,j,k);
                    calc_cell_index(idx2,m2);

                    // add to list of cell neighbors
                    box.cell_neighbors[m1][box.cell_num_neighbors[m1]] = m2;
                    box.do_force_calc[m1][box.cell_num_neighbors[m1]] = box.doNeighborCalc(idx1,idx2);
                    box.cell_num_neighbors[m1]++;

                }
            }
        }
    }

    // now construct the array of list heads and the linked list
    box.head = new int[box.numCell];
    box.bead_list = new int[box.numBeads];

    // and finally, build the lsits
    build_bead_lists();

}


void Sim::build_bead_lists() {
    int i,m;
	
    for (i=0;i<box.numCell;i++) {
        box.head[i] = EMPTY;
    }
 
    for (i=0;i<box.numBeads;i++) {
        box.bead_list[i] = EMPTY;
    }

    for (i=0;i<box.numBeads;i++) {
	
        m = mapToGrid(xs[i]);
        box.bead_list[i] = box.head[m];
        box.head[m] = i;

    }

}


/******************************************************************************************/


/*********************************INTEGRATOR FUNCTIONS*************************************/

void Sim::update_vel_half(NESI::vector *vs, NESI::vector *fs, double dt) {
    
    int i,j;
    double dtHalf = 0.5*dt;

    for (i=0;i<box.numBeads;i++) {
        for (j=0;j<3;j++) {
            vs[i][j] += dtHalf*fs[i][j]*invMasses[types[i]];
        }
    }
}

void Sim::update_pos(NESI::vector *xs, NESI::vector *vs, double dt) {

    int i;
    double randVel;

    // loop through atoms and update velocities
    for (i=0;i<box.numBeads;i++) {

        // update positions
        xs[i][0] += box.dt*vs[i][0];
        xs[i][1] += box.dt*vs[i][1];
        xs[i][2] += box.dt*vs[i][2];

        // handle PBCs
        xs[i][0] = fmod(xs[i][0],box.x);
        if (xs[i][0] < 0.0) xs[i][0] += box.x;
    
        xs[i][1] = fmod(xs[i][1],box.y);
        if (xs[i][1] < 0.0) xs[i][1] += box.y;
    

        // handle walls if necessary
        if (box.hasWallZ) {
            // if a particle has left the box
            if (xs[i][2] >= box.z) {

                // specular reflection
                xs[i][2] = 2.0*box.z - xs[i][2];
                vs[i][2] *= -1.0;

                // draw new velocity parallel to wall, i.e. diffuse reflection
                vs[i][0] = calc_random_gaussian(0.0,box.velVarTop*sqrtInvMasses[types[i]]);
                vs[i][1] = calc_random_gaussian(0.0,box.velVarTop*sqrtInvMasses[types[i]]);
            }
            if (xs[i][2] < 0.0) {
                
                // specular reflection
                xs[i][2] *= -1.0;
                vs[i][2] *= -1.0;

                // draw new velocity
                vs[i][0] = calc_random_gaussian(0.0,box.velVar*sqrtInvMasses[types[i]]);
                vs[i][1] = calc_random_gaussian(0.0,box.velVar*sqrtInvMasses[types[i]]);
            }
        }
        else {
            xs[i][2] = fmod(xs[i][2],box.z);
            if (xs[i][2] < 0.0) xs[i][2] += box.z;
        }

    }
}


// apply canonical sampling velocity rescaling (CSVR) thermostat near the lower wall
// and Andersen thermostat near top wall
void Sim::scale_wall_velocities() {

    int i,j,count;
    double var,kinetic,dK,newK,targetK,scale,height;

    // get kinetic energy and find out how we need to scale
    //#pragma omp parallel for reduction(+:kinetic)
    kinetic = 0;
    count = 0;
    for (i=0;i<box.numBeads;i++){
        // if near lower wall, add to count/kinetic totals for CSVR thermostat
        if (xs[i][2] <= box.thermoHeight) {
            kinetic += masses[types[i]]*vec_dot(vs[i],vs[i]);
            count++;
        }
        // Andersen thermostat at upper wall (nothing to do after this step for upper wall)
        else if (xs[i][2] >= (box.z-box.thermoHeightTop)) {
            if (RanGen->Random() <= box.collisionProb) {
                var = box.velVarTop*sqrtInvMasses[types[i]];
                vs[i][0] = calc_random_gaussian(0.0,var);
                vs[i][1] = calc_random_gaussian(0.0,var);
                vs[i][2] = calc_random_gaussian(0.0,var);    
            }  
        }
    }
    kinetic *= 0.5;

    // get the target kinetic energy for the given number of particles
    targetK = 1.5*count*box.tempSet;

    // take a step in the stochastic process
    dK = (targetK-kinetic)*(box.dt/box.rescaleTime) + 2.0*sqrtl(targetK*kinetic*box.dt/(3.0*count*box.rescaleTime))*calc_random_gaussian(0.0,1.0);
    newK = kinetic + dK;

    // find the scaling parameter
    scale = sqrtl(newK/kinetic);
    
    // scale the velocities
    for (i=0;i<box.numBeads;i++) {
        if (xs[i][2] <= box.thermoHeight) {
            for (j=0;j<3;j++) {
                vs[i][j] *= scale;
            }
        }
    }
    
}

// apply canonical sampling velocity rescaling (CSVR) thermostat
// currently unused
void Sim::scale_velocities() {

    int i,j;
    double kinetic,dK,newK,scale;

    // get kinetic energy and find out how we need to scale
    //#pragma omp parallel for reduction(+:kinetic)
    kinetic = 0;
    for (i=0;i<box.numBeads;i++){
        kinetic += masses[types[i]]*vec_dot(vs[i],vs[i]);
    }
    kinetic *= 0.5;

    // take a step in the stochastic process
    dK = (box.targetK-kinetic)*(box.dt/box.rescaleTime) + 2.0*sqrtl(box.targetK*kinetic*box.dt/(3.0*box.numBeads*box.rescaleTime))*calc_random_gaussian(0.0,1.0);
    newK = kinetic + dK;

    // find the scaling parameter
    scale = sqrtl(newK/kinetic);
    
    // scale the velocities
    for (i=0;i<box.numBeads;i++) {
        for (j=0;j<3;j++) {
            vs[i][j] *= scale;
        }
    }
}


// zero out the momentum near the wall - plays nice with the CSVR thermostat being used
void Sim::zero_wall_momentum() {

    int i,j,count;
    double kin,newKin,scale;
    NESI::vector pTot;

    // get momentum, keeping track of the kinetic energy
    count = 0;
    kin = 0.0;
    pTot[0] = pTot[1] = pTot[2] = 0.0;
    for (i=0;i<box.numBeads;i++){
        if (xs[i][2] < box.pinHeight) {
            pTot[0] += masses[types[i]]*vs[i][0];
            pTot[1] += masses[types[i]]*vs[i][1];
            pTot[1] += masses[types[i]]*vs[i][2];
            kin += masses[types[i]]*vec_dot(vs[i],vs[i]);
            count++;
        }
    }
    
    // calculate kinetic energy and average momentm
    pTot[0] /= count;
    pTot[1] /= count;
    pTot[2] /= count;
    newKin = 0.0;
    
    // shift the velocities, keeping track of the new kinetic energy
    for (i=0;i<box.numBeads;i++) {
        if (xs[i][2] < box.pinHeight) {
            vs[i][0] -= invMasses[types[i]]*pTot[0];
            vs[i][1] -= invMasses[types[i]]*pTot[1];
            vs[i][2] -= invMasses[types[i]]*pTot[2];
            newKin += masses[types[i]]*vec_dot(vs[i],vs[i]);
        }
    }
    
    // scale the velocities back to the same kinetic energy we had before
    scale = sqrtl(kin/newKin);
    for (i=0;i<box.numBeads;i++) {
        if (xs[i][2] < box.pinHeight) {
            vs[i][0] *= scale;
            vs[i][1] *= scale;
            vs[i][2] *= scale;
        }
    }



}

// function to swap particle identities near the walls in order to maintain chemical potential difference at top/bottom
// note that this assumes only two particle types
void Sim::swap_particles() {

    // build vector of atom ids for upper and lower baths
    std::vector<int> upperIds;
    std::vector<int> lowerIds;
    int i,j,k,newType,factor,m1,m2;
    double newPotential,deltaU,boltz;
    
    
    // get the particles that are in the upper/lower baths
    for (i=0;i<box.numBeads;i++) {   
        if (xs[i][2] <= box.thermoHeight) {
            lowerIds.push_back(i);
        }
        else if (xs[i][2] >= (box.z-box.thermoHeightTop)) {
            upperIds.push_back(i);
        }
    }

    
    // first do the upper bath
    // choose a particle at random and figure out what the new type should be
    i = upperIds.at((int) floor(RanGen->Random()*upperIds.size()));
    newPotential = 0.0;
    newType = (types[i]+1)%2;
    factor = types[i]-newType;
    
    // now calculate the new energy using the cell list structure
    // get the grid cell the atom is in
    m1 = mapToGrid(xs[i]);
    j = box.head[m1];
    
    // interactions within the current cell
    while (j != EMPTY) {
        newPotential += calc_position_energy(i,j,newType);
        j = box.bead_list[j];
    }

    for (k=0;k<box.cell_num_neighbors[m1];k++) {

        // get this cell's placement
        m2 = box.cell_neighbors[m1][k];
        
        // loop through the atoms in this cell
        j = box.head[m2];
        while (j != EMPTY) {
            newPotential += calc_position_energy(i,j,newType);
            j = box.bead_list[j];
        }
    }
    
    // get the boltzmann factor
    deltaU = newPotential - 2.0*potentials[i];
    boltz = exp((-1.0*deltaU + box.deltaMuTop*factor)/box.tempSetTop);
    
    // determine if the move is accepted
    if (boltz >= 1.0 || RanGen->Random() < boltz) {
        types[i] = newType;
        box.upperAccept += 1;
    }
    
    
    
    // now do the lower bath
    // choose a particle at random and figure out what the new type is
    i = lowerIds.at((int) floor(RanGen->Random()*lowerIds.size()));
    newPotential = 0.0;
    newType = (types[i]+1)%2;
    factor = types[i]-newType;
    
    // now calculate the new energy using the cell list structure
    // get the grid cell the atom is in
    m1 = mapToGrid(xs[i]);
    j = box.head[m1];
    
    // interactions within the current cell
    while (j != EMPTY) {
        newPotential += calc_position_energy(i,j,newType);
        j = box.bead_list[j];
    }

    for (k=0;k<box.cell_num_neighbors[m1];k++) {

        // get this cell's placement
        m2 = box.cell_neighbors[m1][k];
        
        // loop through the atoms in this cell
        j = box.head[m2];
        while (j != EMPTY) {
            newPotential += calc_position_energy(i,j,newType);
            j = box.bead_list[j];
        }
    }
    
    // get the boltzmann factor
    deltaU = newPotential - 2.0*potentials[i];
    boltz = exp((-1.0*deltaU + box.deltaMuBottom*factor)/box.tempSet);
    
    // determine if the move is accepted
    if (boltz >= 1.0 || RanGen->Random() < boltz) {
        types[i] = newType;
        box.lowerAccept += 1;
    }
    
}




/***************************************************************************************************/


/********************************* DATA AND FORCE CALCULATIONS**************************************/



// calculate slab data (parallelized)
void Sim::calc_slab_data_omp() {

    #pragma omp parallel
    {

        int i,j,k,l,idx,tid,start,stop,chunk,offset,a,b,c,aa,bb,cc,m1,m2,atomIdx,histIdx,slabIdx;
        shape virialTens,kineticTens;
        NESI::vector p,virFlux,enerFlux;
        //int **counts;
       //counts = new int *[box.numTypes];
        //NESI::vector **mtms;
        //mtms = new NESI::vector *[box.numTypes];
        double *pots,*typeEners;
        shape *kins,*virs;
        NESI::vector *fluxes;
        double E,x,y,z,kineticScalar,histPos;
        bool continueCalc;
        std::vector<int *> counts;
        std::vector<NESI::vector *> mtms;

        // get thread id
        #if defined(_OPENMP)
        tid = omp_get_thread_num();
        #else
        tid=0;
        #endif
        

        // get pointers to this thread's area of the array(s)
        for (i=0;i<box.numTypes;i++) {
            counts.push_back(slab_counts[i]+(tid*box.numSlabs));
            mtms.push_back(slab_mtms[i]+(tid*box.numSlabs));
            //counts[i]=slab_counts[i]+(tid*box.numSlabs);
            //mtms[i]=slab_mtms[i]+(tid*box.numSlabs);
        }
        kins=slab_kins+(tid*box.numSlabs);
        pots=slab_pots+(tid*box.numSlabs);
        virs=slab_virs+(tid*box.numSlabs);
        fluxes=slab_heat_fluxes+(tid*box.numSlabs);
        typeEners=typeEnergies+(tid*box.numTypes);

        // zero out the data, but only if it's not the first thread since the
        // first set of values is the accumulator which gets averaged and written later on
        if (tid != 0) {
            for (i=0;i<box.numTypes;i++) {
                memset(counts[i], 0, box.numSlabs*sizeof(*counts[i]));
                memset(mtms[i], 0, box.numSlabs*sizeof(*mtms[i]));
            }
            memset(kins, 0, box.numSlabs*sizeof(*kins));
            memset(pots, 0, box.numSlabs*sizeof(*pots));
            memset(virs, 0, box.numSlabs*sizeof(*virs));
            memset(fluxes, 0, box.numSlabs*sizeof(*fluxes));
        }

        // loop through all atoms and get slab counts, kinetic energy, pressure, velocity, and potential
        for (i=0;i<box.numBeads;i+=nThreads) {
           
            // get the atom id
            j = i+tid;
            if (j >= box.numBeads) break;
            
            // get which slab its in
            idx = floor(xs[j][2]/box.slabWidth);
            
            // calculate counts, kinetic energy, momentum, and virials sum
            counts[types[j]][idx]++;
            vec_direct(kineticTens,vs[j],vs[j]);
            shape_scalar_mult(kineticTens,0.5*masses[types[j]]); 
            shape_add(kins[idx],kineticTens);
            vec_scalar_mult(p,vs[j],masses[types[j]]);
            vec_add(mtms[types[j]][idx],p);
            pots[idx] += potentials[j];

            // virial contribution to stress (factor of 1/2 prevents double-counting)
            shape_mult(virialTens,virials[j],0.5);
            shape_add(virs[idx],virialTens);

            // virial contribution to heat flux
            vec_dot_shape(virFlux,vs[j],virialTens);
            vec_add(fluxes[idx],virFlux);

            // potential and kinetic contribution to heat flux
            kineticScalar = shape_trace(kineticTens);
            vec_scalar_mult(enerFlux,vs[j],kineticScalar+potentials[j]);
            vec_add(fluxes[idx],enerFlux);

        }

        // synch up the threads
        #if defined(_OPENMP)
        #pragma omp barrier
        #endif
        
        // now reduce the arrays
        chunk = 1 + (box.numSlabs/nThreads);  // how many slabs for each thread to reduce
        start = tid * chunk;
        stop = start + chunk;
        if (stop > box.numSlabs) stop = box.numSlabs;

        // each threads now loops through the each copy of the slab data, reducing the specific slabs assigned to it
        // IMPORTANT - skip the first section of the array or you'll double count things!
        // I'm not dealing with normalizing the pressures/velocities/virials/etc. as these are easily done in post-processing
        for (i=1;i<nThreads;i++) {
            
            offset = i*box.numSlabs;
            for (j=start;j<stop;j++) {
                
                for (k=0;k<box.numTypes;k++) {
                    slab_counts[k][j] += slab_counts[k][j+offset];
                    vec_add(slab_mtms[k][j],slab_mtms[k][j+offset]); 
                }
                slab_pots[j] += slab_pots[j+offset];
                shape_add(slab_kins[j],slab_kins[j+offset]);
                shape_add(slab_virs[j],slab_virs[j+offset]);
                vec_add(slab_heat_fluxes[j],slab_heat_fluxes[j+offset]);
            }
        }
           
        // now do the particle insertions since we already have the threads up and running
        // each thread handles a bunch of slabs (important to space them out so that you dont
        // end up with only 1-2 threads handling the entire liquid region)
        for (j=0;j<box.numSlabs;j+=nThreads) {

            // get slab index
            slabIdx = j+tid;
            if (slabIdx >= box.numSlabs) break;

            // for each insertion
            for (i=0;i<box.numInsertionsPerStep;i++) {
                
                //reset energy
                memset(typeEners, 0.0, box.numTypes*sizeof(*typeEners));

                // generate random position in the slab for a trial insertion
                x = box.x*threadRanGen[tid].Random();
                y = box.y*threadRanGen[tid].Random();
                z = box.slabWidth*(threadRanGen[tid].Random()+slabIdx); // make sure to adjust the proper slab height
      
                // find out which cell the position is in
                m1 = mapToGrid(x,y,z);

                // do all intra-cell interactions
                atomIdx = box.head[m1];
                while (atomIdx != EMPTY) {
                    calc_position_energy(x,y,z,atomIdx,typeEners);
                    atomIdx = box.bead_list[atomIdx];
                }

                // see if we're above the energy threshold - if any of the interaction energies is below the quite energy, keep going
                continueCalc = false;
                for (k=0;k<box.numTypes;k++) {
                    if (typeEners[k] <= quitEnergy) {
                        continueCalc = true;
                        break;
                    }
                }

                // loop through atoms in nearby cells
                if (continueCalc) {
                    for (k=0;k<box.cell_num_neighbors[m1];k++) {

                        // get this cell's placement
                        m2 = box.cell_neighbors[m1][k];
                        
                        // loop through the atoms in this cell
                        atomIdx = box.head[m2];
                        while (atomIdx != EMPTY) {
                            
                            // get energy, tracking contributions from each type separately
                            // the epsilon parameters are not included in this function, that comes later
                            calc_position_energy(x,y,z,atomIdx,typeEners);  // only works if sigma's and rcuts are all equal!!!!
                            
                            // update atom
                            atomIdx = box.bead_list[atomIdx];
                        }

                        // see if the energy is high enough and we should exit the loop
                        if (k < box.cell_num_neighbors[m1]-1) {
                            continueCalc = false;
                            for (l=0;l<box.numTypes;l++) {
                                if (typeEners[l] <= quitEnergy) {
                                    continueCalc = true;
                                    break;
                                }
                            }
                        }
                        if (!continueCalc) break;
                    }
                }

                // for each LJ type
                for (k=0;k<box.numTypes;k++) {
                    
                    // find placement in the histogram and record it
                    histPos = (typeEners[k]-minEnergy)/binWidth;
                    if (histPos > numBins) {
                        histIdx = numBins;
                    }
                    else {
                        histIdx = (int) floor(histPos);
                    }
                    slab_widom_hist[k][slabIdx][histIdx]++;
                
                }
            }
        }
    } // end of parallel section
}


// wrapper function to call the various force functions
void Sim::calc_forces(bool doVirials,bool doEnergy) {

    // calculate nonbonded forces
    calc_nonbond_forces_omp(doVirials,doEnergy);
   
    // force to pin the liquid on the lower wall
    if (box.hasWallZ) {
        calc_pin_forces();
        calc_wall_forces(doVirials,doEnergy);
    }

}


// force, virial, and energy calculation - used openmp
void Sim::calc_nonbond_forces_omp(bool doVirials, bool doEnergy) {

    // fire up the threads
    #pragma omp parallel
    {
    
        int tid,i,j,k,l,m1,m2,b1,b2,start,stop,offset;
        NESI::vector *frc;
        NESI::shape *vir;
        double *pot;

        // get thread id 
        #if defined(_OPENMP)
        tid = omp_get_thread_num();
        #else
        tid=0;
        #endif

        // get pointer to start of force array and set it to zero
        frc = fs + (tid*box.numBeads);
        memset(frc, 0.0, box.numBeads*sizeof(*frc));

        vir = virials + (tid*box.numBeads);
        if (doVirials) memset(vir, 0.0, box.numBeads*sizeof(*vir));  // reset virials if necessary

        pot = potentials + (tid*box.numBeads);
        if (doEnergy) memset(pot,0.0,box.numBeads*sizeof(*pot));  // reset energies if necessary

        // loop through cells
        for (i=0;i<box.numCell;i+=nThreads) {

            // get the cell id for this iteration
            m1 = i + tid;
            if (m1 >= box.numCell) break;
            
            // move along if there's no beads in the cell
            if (box.head[m1] == EMPTY) continue;

            // first do intra-cell interactions with a double loop over the cell list
            b1 = box.head[m1];
            while (b1 != EMPTY) {
                b2 = box.head[m1];
                while (b2 != EMPTY) {
                    if (b1 > b2) {
                        calc_pair_force_omp(b1,b2,doVirials,doEnergy,frc,vir,pot);
                    }
                    b2 = box.bead_list[b2];
                }
                b1 = box.bead_list[b1];
            }
            
            // now loop through the cell's neighbors
            for (j=0;j<box.cell_num_neighbors[m1];j++) {
                
                // only do half the neighbors (by Newton's law)
                if (!box.do_force_calc[m1][j]) continue;
                
                // get cell index
                m2 = box.cell_neighbors[m1][j];

                // move along if the cell is empty
                if (box.head[m2] == EMPTY) continue;

                // calculate forces for each pair of beads
                b1 = box.head[m1];
                while (b1 != EMPTY) {
                    b2 = box.head[m2];
                    while (b2 != EMPTY) {
                        calc_pair_force_omp(b1,b2,doVirials,doEnergy,frc,vir,pot);
                        b2 = box.bead_list[b2];
                    }
                    b1 = box.bead_list[b1];
                }
            }
        }
    
        #pragma omp barrier // sync threads
       
        // now combine the forces
        // get equal chunks of atoms to handle
        i = 1 + (box.numBeads/nThreads);
        start = tid*i;
        stop = start + i;
        if (stop > box.numBeads) stop = box.numBeads;
        
        // now reduce to the first spot in the array
        for (i=1;i<nThreads;i++) {
           
            // offset for the next block of data
            offset = i*box.numBeads;

            // for each atom this thread is responsible for
            for (j=start;j<stop;j++) {

                // add the data into the zeroth array
                fs[j][0] += fs[offset+j][0];
                fs[j][1] += fs[offset+j][1];
                fs[j][2] += fs[offset+j][2];

                // reduce virials if necessary
                if (doVirials) {
                    virials[j][0] += virials[offset+j][0];
                    virials[j][1] += virials[offset+j][1];
                    virials[j][2] += virials[offset+j][2];
                    virials[j][3] += virials[offset+j][3];
                    virials[j][4] += virials[offset+j][4];
                    virials[j][5] += virials[offset+j][5];
                }

                if (doEnergy) {
                    potentials[j] += potentials[offset+j];
                }
            }
        }

    } // end of parallel section
}



// force, virial, and energy calculation - used openmp (embarassingly parallel)
void Sim::calc_wall_forces(bool doVirials, bool doEnergy) {

    int tid,i,idx;
    
    // loop through atoms
    for (i=0;i<box.numBeads;i++) {

        // calculate the wall force... all wall energy (and virial contribution) is stored in the atoms
        calc_wall_force(i,doVirials,doEnergy);
    }
}


// force to pin the liquid to the lower wall, ignores virial contributions
void Sim::calc_pin_forces() {

    int i;
    
    // loop through the atoms and apply pinning force if in the specified region
    for (i=0;i<box.numBeads;i++) {
        if (xs[i][2] < box.pinHeight) fs[i][2] -= box.pinningForce;
    }
}



// calculates interaction energy between existing atom (atomIdx) with a new LJ particle with given position (x,y,z) and stores it in the array
// this only works if all particles have same diameter and cutoff
void Sim::calc_position_energy(double x, double y, double z, int atomIdx, double *energies) {
    double sigma,cutoff,distSqr,invR2,invR6,invR12,dist,distMinRc,distMinRcSqr,factor;
    int t,k;
    NESI::vector pos,rnear;

    // get atom types and resulting parameters
    t = types[atomIdx];
    sigma = sigs[t][t];
    cutoff = rcuts[t][t];

    // get distance vector between the atoms
    pos[0] = x;
    pos[1] = y;
    pos[2] = z;
    nearest_image_dist(rnear,xs[atomIdx],pos);
    distSqr = vec_dot(rnear,rnear);

    // ignore if the atoms are beyond the cutoff
    if(cutoff*cutoff < distSqr) {
        for (k=0;k<box.numTypes;k++) {
            energies[k] += 0.0;
        }
    }
    // if we're above the spline cutoff
    else if (rs*rs < distSqr) {
        dist = sqrtl(distSqr);
        distMinRc = dist - cutoff;
        distMinRcSqr = distMinRc*distMinRc;     
        factor = xi1*distMinRc*distMinRcSqr + xi2*distMinRcSqr;
        for (k=0;k<box.numTypes;k++) {
            energies[k] += eps[t][k]*factor;
        }
    }
    // otherwise do the full LJ potential
    else {
        invR2 = sigma*sigma/distSqr;
        invR6 = invR2*invR2*invR2;
        invR12 = invR6*invR6;
        factor= 4.0*(invR12-invR6);
        for (k=0;k<box.numTypes;k++) {
            energies[k] += eps[t][k]*factor;
        }
    }

}

// calculates interaction energy between existing atom (atomIdx) with a new LJ particle with given position (x,y,z) and stores it in the array
// epsilons are not considered and will be dealt with outside this function
double Sim::calc_position_energy(int i, int j, int t1) {
    
    if (i == j) {
        return 0.0;
    }
    
    double sigma,cutoff,epsilon,distSqr,invR2,invR6,invR12,dist,distMinRc,distMinRcSqr;
    int t2;
    NESI::vector rnear;

    // get atom types and resulting parameters
    t2 = types[j];
    sigma = sigs[t1][t2];
    cutoff = rcuts[t1][t2];
    epsilon = eps[t1][t2];

    // get distance vector between the atoms
    nearest_image_dist(rnear,xs[i],xs[j]);
    distSqr = vec_dot(rnear,rnear);

    // ignore if the atoms are beyond the cutoff
    if(cutoff*cutoff < distSqr) {
        return 0.0;
    }
    // if we're above the spline cutoff
    else if (rs*rs < distSqr) {
        dist = sqrtl(distSqr);
        distMinRc = dist - cutoff;
        distMinRcSqr = distMinRc*distMinRc;
        return epsilon*(xi1*distMinRc*distMinRcSqr + xi2*distMinRcSqr);
    }
    // otherwise do the full LJ potential
    else {
        invR2 = sigma*sigma/distSqr;
        invR6 = invR2*invR2*invR2;
        invR12 = invR6*invR6;
        return 4.0*epsilon*(invR12-invR6);
    }

}



// calculates the force on atom i due to atom j and stores it in the given pointer array
void Sim::calc_pair_force_omp(int i, int j, bool doVirials, bool doEnergy, NESI::vector *frc, shape *vir, double *pot) {
    double halfE,factor,epsilon,sigma,cutoff,distSqr,dist,distMinRc,invR1,invR2,invR6,invR8,invR12,invR14;
    int t1,t2;
    NESI::vector dr,force;
    shape virial;

    // get atom types and resulting parameters
    t1 = types[i];
    t2 = types[j];
    epsilon = eps[t1][t2];
    sigma = sigs[t1][t2];
    cutoff = rcuts[t1][t2];

    // get distance vector between the atoms and store in the force vector (they're proportional!)
    nearest_image_dist(dr,xs[i],xs[j]); // dr = r_i - r_j
    distSqr = vec_dot(dr,dr);

    // ignore if the atoms are beyond the cutoff
    if(cutoff*cutoff < distSqr) {
        return;
    }
    
    // if we're above the spline cutoff
    if (rs*rs < distSqr) {
        dist = sqrtl(distSqr);
        distMinRc = dist - cutoff;
        factor = -1.0*epsilon*(3.0*xi1*distMinRc*distMinRc + 2.0*xi2*distMinRc)/dist;
        if (doEnergy) halfE = 0.5*epsilon*(xi1*distMinRc*distMinRc*distMinRc + xi2*distMinRc*distMinRc);
    }
    else {
        // calculate the factors we need
        invR2 = sigma*sigma/distSqr;
        invR6 = invR2*invR2*invR2;
        invR12 = invR6*invR6;
        invR8 = invR6*invR2;
        invR14 = invR12*invR2;
        invR1 = sqrtl(invR2);
        factor = 24.0*epsilon*(2.0*invR14-invR8);
        if (doEnergy) halfE = 2.0*epsilon*(invR12-invR6); // each particle gets half the energy
    } 
    
    // get the force    
    force[0] = dr[0]*factor;
    force[1] = dr[1]*factor;
    force[2] = dr[2]*factor;

    // update the forces
    vec_add(frc[i],force);
    vec_subtr(frc[j],force);

    // do virials if necessary
    if (doVirials) {
        vec_direct(virial,force,dr);
        shape_add(vir[i],virial);
        shape_add(vir[j],virial);
    }

    //do energy if necessary
    if (doEnergy) {
        pot[i] += halfE;
        pot[j] += halfE;
    }

}


// calculates the force on atom i due to atom j and store in the given pointer array
void Sim::calc_wall_force(int i, bool doVirials, bool doEnergy) {
    double halfE,force,epsilon,sigma,dist,distMinRc,invR1,invR2,invR4,invR5,invR10,invR11;
    int t1;


    // get distance between wall and atom
    dist = xs[i][2];

    // if we're above the cutoff, just continue
    if (dist > box.wallCut) {
        return;
    }

    // otherwise, get atom types and resulting parameters
    t1 = types[i];
    sigma = sigs[t1][t1];
   

    // if we're above the spline cutoff
    if (rs < dist) {
        distMinRc = dist - box.wallCut;
        force = -1.0*box.wallEps*(3.0*xi1Wall*distMinRc*distMinRc + 2.0*xi2Wall*distMinRc);
        if (doEnergy) halfE = 0.5*box.wallEps*(xi1Wall*distMinRc*distMinRc*distMinRc + xi2Wall*distMinRc*distMinRc);
    }
    else {
        // calculate the factors we need
        invR1 = sigma/dist;
        invR2 = invR1*invR1;
        invR4 = invR2*invR2;
        invR5 = invR1*invR4;
        invR10 = invR5*invR5;
        invR11 = invR10*invR1;
        force = 6.4*PI*box.wallEps*(invR11/3.0-invR5);
        if (doEnergy) halfE = 0.8*PI*box.wallEps*((2.0/15.0)*invR10-invR4); // each particle gets half the energy - the wall gets the other half (I guess!)
    } 
    
    // update the forces
    fs[i][2] += force;

    // do virials if necessary
    if (doVirials) {
        //virials[i][2] += force;
    }

    //do energy if necessary
    if (doEnergy) {
        //potentials[i] += halfE;
    }

}

void Sim::nearest_image_dist(NESI::vector &r, NESI::vector r1, NESI::vector r2) { //r1 - r2 using the nearest image convention
    
    double dx,dy,dz;
    double xhalf = box.x/2;
    double yhalf = box.y/2;
    double zhalf = box.z/2;

    dx = r1[0] - r2[0];
    dy = r1[1] - r2[1];
    dz = r1[2] - r2[2];

    dx = dx > xhalf ? dx - box.x : dx; 
    dx = dx < -xhalf ? dx + box.x : dx; 
    dy = dy > yhalf ? dy - box.y : dy; 
    dy = dy < -yhalf ? dy + box.y : dy; 
    dz = dz > zhalf ? dz - box.z : dz; 
    dz = dz < -zhalf ? dz + box.z : dz; 

    r[0] = dx; r[1] = dy; r[2] = dz;

}



/**************************************************************************************************************/


/**********************************DATA OUTPUT AND VISUALIZATION***********************************************/


// write checkpoint file and widom histograms
void Sim::print_configuration() {
    int i,j,k;
    char tt[2001];
    int start,end,anchor,rigid,issite,hassites;
    FILE *atom_file = fopen("beads.dat","w"); //Location of all atom files 
    fprintf(atom_file, "%d\n",box.numBeads);
    fprintf(atom_file, "%d\n",box.numTypes);
    for (i=0;i<box.numBeads;i++) {
        fprintf(atom_file,"%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n",types[i],xs[i][0],xs[i][1],xs[i][2],vs[i][0],vs[i][1],vs[i][2]);
    }

    fclose(atom_file);

    // write the widom histograms
    for (i=0;i<box.numTypes;i++) {
        char filename[32];
        sprintf(filename,"slab_widom_hist_%d.txt",i); 
        FILE *hist = fopen(filename,"a");
        for (j=0;j<box.numSlabs;j++) {
            fprintf(hist,"%d",j);
            for (k=0;k<numBins+1;k++){
                fprintf(hist,"\t%ld",slab_widom_hist[i][j][k]);
            }
            fprintf(hist,"\n");
        }
        fclose(hist);
    }

    // reset the histograms
    for (i=0;i<box.numTypes;i++) {
        for (j=0;j<box.numSlabs;j++) {
            memset(slab_widom_hist[i][j], 0, (numBins+1)*sizeof(*slab_widom_hist[i][j]));
        }
    }

}

// handle boundary conditions
void Sim::wrap_atoms() {
    
    int i;

    for (i=0;i<box.numBeads;i++) {
            
        xs[i][0] = fmod(xs[i][0],box.x);
        if (xs[i][0] < 0.0) xs[i][0] += box.x;
    
        xs[i][1] = fmod(xs[i][1],box.y);
        if (xs[i][1] < 0.0) xs[i][1] += box.y;
    
        xs[i][2] = fmod(xs[i][2],box.z);
        if (xs[i][2] < 0.0) xs[i][2] += box.z;
    
    }
}

void Sim::printXYZ() {
    //Print all the molecules then print out rings to visualize disk-like "atoms"
    
    int i,j;
    FILE *dump;
    dump = fopen("config.xyz","a");
    char typeChars[] = "CNOBF"; // convert type integers to atom abbreviations for visualization

    fprintf(dump,"%d\nConfiguration\n",box.numBeads);
    for (i=0;i<box.numBeads;i++) {
        fprintf(dump,"%c\t%f\t%f\t%f\n",typeChars[types[i]],xs[i][0],xs[i][1],xs[i][2]);
    }
    fclose(dump);
}



// write slab data, averaging as we go, then reset the arrays
void Sim::write_slab_data() {
   
    int i,j;
    double nFrames = floor((double) box.writeInterval/box.calcInterval);
    //calc_slab_data_omp();
    double slabVol = box.x*box.y*box.slabWidth;

  
    // virial tensor
    FILE *vir = fopen("slab_virs.txt","a");
    fprintf(vir,"%d",box.step);
    for (i=0;i<box.numSlabs;i++){
        shape_scalar_mult(slab_virs[i],1.0/(nFrames*slabVol));
        fprintf(vir,"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",slab_virs[i][0],slab_virs[i][1],slab_virs[i][2],slab_virs[i][3],slab_virs[i][4],slab_virs[i][5]); // full virial tensor
    }
    fprintf(vir,"\n");
    fclose(vir);


  
    // kinetic energy density (in tensor form)
    FILE *kin = fopen("slab_kins.txt","a");
    fprintf(kin,"%d",box.step);
    for (i=0;i<box.numSlabs;i++){
        shape_scalar_mult(slab_kins[i],1.0/(nFrames*slabVol));
        fprintf(kin,"\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf",slab_kins[i][0],slab_kins[i][1],slab_kins[i][2],slab_kins[i][3],slab_kins[i][4],slab_kins[i][5]); // full virial tensor
    }
    fprintf(kin,"\n");
    fclose(kin);

    // species mass densities
    for (j=0;j<box.numTypes;j++) {
        char filename[32];
        sprintf(filename,"slab_dens_%d.txt",j); 
        FILE *dens = fopen(filename,"a");
   
        fprintf(dens,"%d",box.step);
        for (i=0;i<box.numSlabs;i++){
            fprintf(dens,"\t%lf",(slab_counts[j][i]*masses[types[j]])/(nFrames*slabVol));
        }
        fprintf(dens,"\n");
        fclose(dens);
    }

    // potential energy densities
    FILE *pot = fopen("slab_pots.txt","a");
    fprintf(pot,"%d",box.step);
    for (i=0;i<box.numSlabs;i++){
        fprintf(pot,"\t%lf",slab_pots[i]/(nFrames*slabVol));
    }
    fprintf(pot,"\n");
    fclose(pot);

   
    // species momentum densities
    for (j=0;j<box.numTypes;j++) {
        char filename[32];
        sprintf(filename,"slab_mtms_%d.txt",j); 
        FILE *mtms = fopen(filename,"a");
   
        fprintf(mtms,"%d",box.step);
        for (i=0;i<box.numSlabs;i++){
            vec_scalar_mult(slab_mtms[j][i],1.0/(nFrames*slabVol));
            fprintf(mtms,"\t%lf\t%lf\t%lf",slab_mtms[j][i][0],slab_mtms[j][i][1],slab_mtms[j][i][2]);
        }
        fprintf(mtms,"\n");
        fclose(mtms);
    }



    // heat fluxes
    FILE *flux = fopen("slab_heat_fluxes.txt","a");
    fprintf(flux,"%d",box.step);
    for (i=0;i<box.numSlabs;i++){
        vec_scalar_mult(slab_heat_fluxes[i],1.0/(nFrames*slabVol));
        fprintf(flux,"\t%lf\t%lf\t%lf",slab_heat_fluxes[i][0],slab_heat_fluxes[i][1],slab_heat_fluxes[i][2]);
    }
    fprintf(flux,"\n");
    fclose(flux);

    // now reset the arrays - only need to worry about the first chunk since the rest are automatically reset
    // for each slab calculation
    for (i=0;i<box.numTypes;i++) {
        memset(slab_counts[i], 0, box.numSlabs*sizeof(*slab_counts[i]));
        memset(slab_mtms[i], 0, box.numSlabs*sizeof(*slab_mtms[i]));
        
    }
    memset(slab_kins, 0, box.numSlabs*sizeof(*slab_kins));
    memset(slab_pots, 0, box.numSlabs*sizeof(*slab_pots));
    memset(slab_virs, 0, box.numSlabs*sizeof(*slab_virs));
    memset(slab_heat_fluxes, 0, box.numSlabs*sizeof(*slab_heat_fluxes));

}

void Sim::dumpData() {
    FILE *out = fopen("data.txt","w");
    int i,j,k;
    fprintf(out,"%d\n%lf\t%lf\t%lf\n",box.numBeads,box.x,box.y,box.z);
    for (i=0;i<box.numBeads;i++) {
        fprintf(out,"%d\t", types[i]);
        fprintf(out,"%lf\t%lf\t%lf\t", xs[i][0], xs[i][1], xs[i][2]);
        fprintf(out,"%lf\t%lf\t%lf\n", vs[i][0], vs[i][1], vs[i][2]);

    }
    fclose(out);

	/*
    for (i=0;i<box.numTypes;i++) {
        char filename[32];
        sprintf(filename,"slab_widom_hist_%d.txt",i); 
        FILE *hist = fopen(filename,"a");
        for (j=0;j<box.numSlabs;j++) {
            fprintf(hist,"%d",j);
            for (k=0;k<numBins+1;k++){
                fprintf(hist,"\t%ld",slab_widom_hist[i][j][k]);
            }
            fprintf(hist,"\n");
        }
        fclose(hist);
    }*/

}


/**************************************************************************************************************/



/******************************************* RNG WRAPPERS *************************************************/


void Sim::calc_random_vector(NESI::vector &b) {
    double R1,R2,R3;
    do { //Generate random unit vector and make it u.
        R1 = (2*RanGen->Random()-1);
        R2 = (2*RanGen->Random()-1);
        R3 = R1*R1+R2*R2;
    } while (R3>=1);
               
    b[0] = 2*sqrtl(1.0-R3)*R1;
    b[1] = 2*sqrtl(1.0-R3)*R2;
    b[2] = 1-2.0*R3;


}


double Sim::calc_random_gaussian(double mu, double sigma) {
    double R1,R2,R3;
    double z;
    do { //Generate random unit vector and make it u.
        R1 = (2*RanGen->Random()-1);
        R2 = (2*RanGen->Random()-1);
        R3 = R1*R1+R2*R2;
    } while (R3>1);
    
    z = sqrtl(-2.0*log(R3)/R3) * R1;
    return z * sigma + mu;

}


/**************************************************************************************************************/

/**************************************************************************************************************/
