/* Serial Code for N Body Simulation */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include "mpi.h"

/*Global Variables*/
const int numiterations; /*number of resultant iterations */
const int simulation_time = 10; /*Time to run simulation for */
const float dt = 0.01; /*time step every iteration*/
const float G = 1;

/*Position Vector*/
typedef struct {
    double x, y, z; /*three dimensional positions */
    double vx,vy,vz; /*three dimensional velocity*/
    double ax,ay,az; /*three dimensional acceleration*/
    double mass;
} Particle;

/*Particle Vector*/
Particle*  particles;   /* Particle array storing structs for each particle, number of particles will be defined in the main function */
int npart;

/*Initialize particles positions, velocity and acceleration vectors */
void initParticles(int num_particles){
    int i;
    particles = malloc(sizeof(Particle)*num_particles);
    for (i=0; i<num_particles; i++) {
	particles[i].x	  = drand48();
	particles[i].y	  = drand48();
	particles[i].z	  = drand48();
	particles[i].mass = 1.0;
	particles[i].vx	  = 0;
	particles[i].vy	  = 0;
	particles[i].vz	  = 0;
	particles[i].ax	  = 0;
	particles[i].ay	  = 0;
	particles[i].az	  = 0;
    }
}

/*Update position, velocity, acceleration on each particle*/
void updateParticles(){
    double t = 0;
    int i = 0;
    int k = 0;
    for (t = 0; t < simulation_time; t+=dt){
	for (i = 0; i < npart; i++){
	    double r2 = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	    for (k = 0; k < 3; k++){
		a[k] = - r[k] / (r2 * sqrt(r2));
	    }
	    for (k = 0; k < 3; k++){
		r[k] += v[k] * dt;
		v[k] += a[k] * dt;
	    }
	}
    }
}

int main(int argc, char *argv[]) {	

    int num_particles = 0;
    num_particles = atoi(argv[1]); /*command line input expected which provides the number of particles */
    npart = num_particles;
    initParticles(num_particles);
    updateParticles();
}
