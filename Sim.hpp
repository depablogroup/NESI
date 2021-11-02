#include "Objects.hpp"
#include "Box.hpp"
#include "mersenne.h"
#include "math_vector.h"
#include <vector>

using namespace NESI;

class Sim {
	public:
        Box box; // box (includes grid and linked list)
        int *types; // atom types
        vector *xs; // positions
        vector *vs; // velocities
        vector *fs; // forces
        shape *virials; // per-particle virials
        double *potentials; // per-particle potentials
        int nThreads;

        // slab data
        int **slab_counts;
        shape *slab_kins;
        double *slab_pots;
        vector **slab_mtms;
        shape *slab_virs;
        vector *slab_heat_fluxes;
        double minEnergy;
        double maxEnergy;
        double quitEnergy;
        double binWidth;
        int numBins;
        double *typeEnergies;
        long int ***slab_widom_hist;
        double minEps;

        CRandomMersenne *RanGen;
        CRandomMersenne *threadRanGen;
        void random_thread(int);

        double **eps; // matrix of lennard-jones epsilon parameters
        double **sigs; // matrix of lennard-jones sigma parameters
        double **rcuts; // matrix of lennard-jones cutoffs
        double *masses; // list of atom type masses
        double *invMasses;  // list of inverse masses (for applying forces more quickly)
        double *sqrtInvMasses;  // square roots of inverse masses

        bool calcVirials;
        bool calcEnergy;

        Sim();
        void MDSim();
        void updateT();
        void openFiles();
        void read_input();
        void initialize_data();
        void initialize_system();
        void init_configuration();
        void init_config_liq();
        void set_velocities();
        void allocate_memory();
        void calc_random_vector(vector &b);
        double calc_random_gaussian(double,double);
        void calc_random_normal_vector(vector &b, vector a);
        void calc_cross_vector(vector &c, vector a, vector b);
        double calc_dist(vector r1, vector r2);
        void nearest_image_dist(vector &r, vector r1, vector r2);
        void PBC_shift(vector &r_shift, vector r);
        void proj_vector(vector &a, vector n);
        void moveGrid();
        void initialize_grid();
        void build_bead_lists();
        void build_bead_lists_alt();

        void calc_slab_data();
        void calc_slab_data_omp();
        void calc_slab_insertions(double *tempArray,int type);

        int mapToGrid(NESI::vector r);
        int mapToGrid(double,double,double);
        int calc_cell_placement(int,int,int);
        void calc_cell_index(lvector &index, int);
        void addBeadToGrid(int);
        void removeBeadFromGrid(int);

        void read_configuration();
        void read_topology();
        void print_configuration();
        void generate_lists();
        void generate_matrices();
        void calc_forces(bool,bool);
        void calc_nonbond_forces(bool);
        void calc_nonbond_forces_omp(bool,bool);
        void calc_nonbond_forces_omp_alt(bool,bool);
        void calc_pin_forces();
        void calc_wall_forces(bool,bool);
        void calc_langevin_forces();
        void update_vel_half(vector *vs, vector *fs, double dt);
        void update_pos(vector *xs, vector *vs, double dt);
        double calc_total_energy();
        double calc_pressure(double);
        double calc_total_nonbond_energy();
        void calc_slab_nonbond_energy();
        double calc_total_kinetic_energy(vector *vs);
        void scale_velocities();
        void scale_wall_velocities();
        void swap_particles();
        void zero_wall_momentum();
        double calc_pair_energy(int,int);
        void calc_position_energy(double,double,double,int,double*);
        double calc_position_energy(int,int,int);
        void calc_pair_force(int,int,bool);
        void calc_wall_force(int,bool,bool);
        void calc_pair_force_omp(int,int,bool,bool, vector *frc, shape *vir, double *pot);
        double calc_nonbond_bead_energy(int);


        void wrap_atoms();
        void shiftCOM();
        void printXYZ();
        void writeEnergy();
        void write_slab_data();
        void dumpData();
        void dumpGrid();

        void printCheckpoint();
        void readCheckpoint();
};
