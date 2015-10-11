/* Serial Code for N Body Simulation */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
//#include "mpi.h"

/*Global Variables*/
const int numiterations; /*number of resultant iterations */
const int simulation_time = 10; /*Time to run simulation for */
const double dt = 0.01; /*time step every iteration*/
const double G = 1;

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
double drand(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}
/*Initialize particles positions, velocity and acceleration vectors */
void initParticles(int num_particles){
    particles = malloc(sizeof(Particle)*num_particles);
    for (int i=0; i<num_particles; i++) {
	particles[i].x	  = drand(0, 1);
	particles[i].y	  = drand(0, 1);
	particles[i].z	  = drand(0, 1);
	particles[i].mass = 1.0;
	particles[i].vx	  = 0;
	particles[i].vy	  = 0;
	particles[i].vz	  = 0;
	particles[i].ax	  = 0;
	particles[i].ay	  = 0;
	particles[i].az	  = 0;
    }
}

/*Compute accelerations experienced by each particle*/
void computeAccelerations(){
    /*calculate the acceleration experienced by the particle due to all other particles in this time step*/
    double r2, r3, acx, acy, acz;
    for (int i = 0; i < npart; i++){
	for(int j = 0; j < npart; j++){
	    if (j != i)
	    {
		/*for every pair (i,j) except (i,i) do : */
		double rx = particles[j].x - particles[i].x; /*relative distance between particle j and i based on x dimension*/
		double ry = particles[j].y - particles[i].y; /*relative distance between particle j and i based on y dimension*/
		double rz = particles[j].z - particles[i].z; /*relative distance between particle j and i based on z dimension*/
		r2 = rx*rx + ry*ry + rz*rz; /*magnitude of the distance squared between the particles i and j using euclidean distance*/
		r3 = r2 * sqrt(r2); /*distance to the power of 3*/
		/*acceleration based on newton's gravitational force equation, masses are unity, and G is 1, could change to define the general equation*/
		acx += rx / r3;
		acy += ry/ r3;
		acz += rz / r3;
	    }
	}
	/*update acceleration values*/
	particles[i].ax = acx;
	particles[i].ay = acy;
	particles[i].az = acz;
	/*For each particle print acceleration in x,y,z dimension folllowed by velocity in x,y,z dimension and lastly the position in x,y,z*/
	printf("Particle %d:", i);
	printf(" %f,%f,%f  %f,%f,%f  %f,%f,%f \n", acx, acy, acz, particles[i].vx, particles[i].vy, particles[i].vz,
		particles[i].x, particles[i].y,particles[i].z);
    }
    printf("\n");
}
/*Update position and velocity for next time step based on acceleration computed earlier in computeAccelerations() for this timestep*/
void updateParticles(){
    for (int i = 0; i < npart; i++){
	/*for all particles partcipating in the simulation */
	/*update velocity vectors in each dimension*/
	/*equation for update derived from book which temp shared and used by several other reference codes*/
	particles[i].vx = particles[i].vx + particles[i].ax * dt;
	particles[i].vy = particles[i].vy + particles[i].ay * dt;
	particles[i].vz = particles[i].vz + particles[i].az * dt;
	/*update position vectors in each dimension*/
	/*equation for update dervied from book which temp shared and used by several other reference codes*/
	particles[i].x = particles[i].x + particles[i].vx * dt;
	particles[i].y = particles[i].y + particles[i].vy * dt;
	particles[i].z = particles[i].z + particles[i].vz * dt;
    }
}

void begin_simulation(){
    for(double t = 0; t < simulation_time; t+=dt){
	/*For every iteration in the simulation, which is controlled by dt, the time step and the end time which is simulation_time*/
	printf("Time step t:\n");
	computeAccelerations();
	updateParticles();
    }
}
int main(int argc, char *argv[]){
    int num_particles = 0;
    num_particles = atoi(argv[1]); /*command line input expected which provides the number of particles */
    npart = num_particles;
    /*initialize the position of particles in space and also their velocity and acceleration vector values*/
    initParticles(num_particles);
    /*begin the simulation with the initialization done previously and using simulation parameters defined in the global constants*/
    begin_simulation();
}
