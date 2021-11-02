#include "math_vector.h"
#include <vector>

using namespace NESI;


class QTensor { //Symmetric, traceless 3x3 quadrupole tensor
    public:
        double val[3][3]; //3x3 matrix of values
        double invdens;   //Inverse density - used for incorporating vectors
        void add_vec(vector v); //Adds the effect of a vector into the matrix
        void remove_vec(vector v); //Removes the effect of a vector from the matrix
        void clear();           //Clears the matrix
        double trQ2();          //Calculates the trace of Q^2
        double trQ3();          //Calculates the trace of Q^3
        QTensor();
};

QTensor::QTensor() {
    int i,j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            val[i][j] = 0.0;
        }
    }
}

void QTensor::add_vec(vector v) {
    int i,j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            val[i][j] += invdens*v[i]*v[j];
        }
        val[i][i] -= invdens/3.0;
    }
}

void QTensor::remove_vec(vector v) {
    int i,j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            val[i][j] -= invdens*v[i]*v[j];
        }
        val[i][i] += invdens/3.0;
    }
   
}

double QTensor::trQ2() {
    int i,j;
    double sum = 0.0;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            sum += val[i][j]*val[j][i];
        }
    }
    return sum;
}

double QTensor::trQ3() {
    int i,j,k;
    double sum = 0.0;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            for (k=0;k<3;k++) {
                sum += val[i][j]*val[j][k]*val[k][i]; 
            }
        }
    }
    return sum;
}

void QTensor::clear() {
    int i,j;
    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            val[i][j] = 0.0;
        }
    }
}

class Bead {
    public:
        vector r; //Position within simulation box
        vector rn; //Unwrapped position
        vector u; //Orientation unit vector
        vector f; //Orientation unit vector
        vector v; //Orientation unit vector
        double vol; //Effective volume of the molecule - will be used in calculating volume fractions, etc
        int type; //Type associated with polymer bead
        int grid; //Grid value associated with polymer bead
        bool isAnchor; //Bead is anchored to a surface of some sort and should not move througout the simulation
        int id;     //Index of bead in array
        bool isRod; //If true, then the bead contributes to the Q tensor. If false, then it does not -- useful for speeding up the simulation if nematic interactions are not necessary
        bool isSite; //If true, then this particle will not be allowed to move independently of its parent 
        bool hasSites; //If true, then the object has sites that must be updated simultaneously when performing Monte Carlo moves.
        Bead *next,*prev; //Pointers to the next and previous value in the linked list
};

class linked_list {

    public:
        Bead *head;
        
        linked_list() {
            head = NULL;
        }

        void pushNode(Bead *head_ref) { //Push node to be the head of the list
            head_ref->next = NULL;
            head_ref->prev = NULL;
            
            head_ref->next = head;
            if(head != NULL) head->prev = head_ref;
            head = head_ref;
        }

        /*
        void insertAfter(Bead *new_node, Bead *old_node) { //Insert new_node in postion after old_node
            Bead *next_node = old_node->next;
            if(old_node != NULL) new_node->next = old_node->next;
            new_node->prev = old_node;
            if(old_node != NULL) old_node->next = new_node;
            if(next_node != NULL) next_node->prev = new_node;
        }
        */

        void deleteNode(Bead *node) { //Remove the Bead from the linked list
            
            Bead *prev_node, *next_node;
            prev_node = node->prev;
            next_node = node->next;

            if (prev_node != NULL) prev_node->next = next_node;
            if (next_node != NULL) next_node->prev = prev_node;
            if (node == head) {
                head = next_node;
            }

            node->next = NULL;
            node->prev = NULL;

        }

        int countLength() {
            int i=0;
            Bead *test;
            test = head;
            while(test!=NULL) {
                i++;
                test = test->next;
            }
            return i;
            
        }
};

class Grid {
    public:    
        //QTensor Q; //Local Q-Tensor order parameter associated with the grid site
        linked_list bead_list; //Linked list of bead indices in this grid
        int x,y,z; //Integers specifying the location of the grid cell in 3D space
        int id;       //Location within the grid array

};

class Bond { //Gaussian bonds - minimum energy is always zero
    public:
        double k_G; //Describes harmonic Gaussian bond potential
        double r0; //Describes minimum energy distance
        int id1;    //This id needs to point to the second one and is how u is defined 
        int id2;    //This id is pointed to 

};

class Angle { //Three body angle potential - Minimum always for 180 degrees
    public:
        double k_b; //Bending energy coefficient
        double theta0;
        int id0; //Center body over which angle is defined
        int id1; //Outer body
        int id2; //Outer body

};

class Twist { //Potential related to a dihedral angle between two disks - not a standard dihedral potential
    public:
        //Vector between atom 2 and 3 defines the main axis for the calculation
        //Vectors between 1 and 2 and 3 and 4 are the normal vectors for the plane
        //If id1=id2 or id3=id4, then we will use the f vector as the normal vector (it is disklike)
        int id1;
        int id2;
        int id3; 
        int id4;
        double a1; //Coefficients for cosine part of the Fourier series
        double a2;
        double a3;
        double a4;

};

class Attachment {  //Rigid bond without any energy -- Defines sites on non-pointlike atoms that other atoms can attach to
    public:
        int id_core; //ID of the core "atom" with sites
        int id_site; //ID of the dummy "atom" acting as the attached site

        //Location of the site relative to the core atom's ufv coordinate frame
        double u_comp; 
        double f_comp;
        double v_comp;
};


