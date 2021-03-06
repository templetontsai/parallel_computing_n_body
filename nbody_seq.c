/* Serial Code for N Body Simulation */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include <time.h>
//#include <mpi.h>

/*Global Variables*/
const double dt = 1; /*time step every iteration*/
const double G = 6.673e-11;


/*Position Vector*/
typedef struct {
    double x, y; /*two dimensional positions */
    double vx,vy; /*two dimensional velocity*/
    double ax,ay; /*two dimensional acceleration*/
    double mass;
} Particle;

/*Particle Vector*/
Particle*  particles;   /* Particle array storing structs for each particle, number of particles will be defined in the main function */
int npart = 0;/*NUmber of Particles*/
int nstep = 0;/*Number of Steps*/


double drand(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

float timediff(struct timespec t1, struct timespec t2)
{  
    if (t1.tv_nsec > t2.tv_nsec) {

        t2.tv_sec -= 1;
        t2.tv_nsec += 1000000000;

    }

    return t2.tv_sec-t1.tv_sec + 0.000000001 * (t2.tv_nsec-t1.tv_nsec);
}


/*Initialize particles positions, velocity and acceleration vectors */
void initParticles(int num_particles){
    particles = malloc(sizeof(Particle)*num_particles);
    int i = 0;
    for (i=0; i<num_particles; i++) {
	particles[i].x	  = drand(0, 1);
	particles[i].y	  = drand(0, 1);
	particles[i].mass = 1.0;
#if 0
	particles[i].vx	  = 0;
	particles[i].vy	  = 0;
	particles[i].vz	  = 0;
#endif
	particles[i].vx	  = drand(0, 1);
	particles[i].vy	  = drand(0, 1);
	particles[i].ax	  = drand(0, 1);
	particles[i].ay	  = drand(0, 1);
    }
}

/*Compute accelerations experienced by each particle*/
void computeAccelerations(){
    /*calculate the acceleration experienced by the particle due to all other particles in this time step*/
    double r2, r3, acx, acy, acz, dv_x, dv_y;
    int i, j;
    for (i = 0; i < npart; i++){
	acx = 0.0;
	acy = 0.0;
	for(j = 0; j < npart; j++){
	    if (j == i)
		continue;

	    /*for every pair (i,j) except (i,i) do : */
	    double rx = particles[j].x - particles[i].x; /*relative distance between particle j and i based on x dimension*/
	    double ry = particles[j].y - particles[i].y; /*relative distance between particle j and i based on y dimension*/
	    r2 = rx*rx + ry*ry; /*magnitude of the distance squared between the particles i and j using euclidean distance*/
	    r3 = r2 * sqrt(r2); /*distance to the power of 3*/
	    //printf("r3 = %g", r3);
	    /*acceleration based on newton's gravitational force equation, masses are unity, and G is 1, could change to define the general equation*/
	    acx += (G * particles[i].mass * particles[j].mass * rx) / r3;
	    acy += (G * particles[i].mass * particles[j].mass * ry)/ r3;


	}
	/*update acceleration values*/
	particles[i].ax = acx;
	particles[i].ay = acy;
	/*For each particle print acceleration in x,y,z dimension folllowed by velocity in x,y,z dimension and lastly the position in x,y,z*/
	//printf("Particle %d:", i);
	//printf(" %g,%g,%g  %g,%g,%g  %g,%g,%g \n", acx, acy, acz, particles[i].vx, particles[i].vy, particles[i].vz,
	//	particles[i].x, particles[i].y,particles[i].z);
	/*for all particles partcipating in the simulation */
	/*update velocity vectors in each dimension*/
	/*equation for update derived from book which temp shared and used by several other reference codes*/
	dv_x = (particles[i].ax * dt)/particles[i].mass;
	dv_y = (particles[i].ay * dt)/particles[i].mass;
	particles[i].vx = particles[i].vx + dv_x;
	particles[i].vy = particles[i].vy + dv_y;
	/*update position vectors in each dimension*/
	/*equation for update dervied from book which temp shared and used by several other reference codes*/
	particles[i].x = particles[i].x + (particles[i].vx + dv_x/2) * dt;
	particles[i].y = particles[i].y + (particles[i].vy + dv_y/2) * dt;
	//printf(" vx = %g, vy = %g, vz = %g \n", particles[i].vx, particles[i].vy, particles[i].vz);
    }
    //printf("\n");
}

void begin_simulation(){
    double t = 0;
    for(t = 0; t < nstep; t++){
	/*For every iteration in the simulation, which is controlled by dt, the time step and the end time which is simulation_time*/
	computeAccelerations();

    }
}
int main(int argc, char *argv[]){
    struct timespec bgn,nd;

    npart = atoi(argv[1]); /*command line input expected which provides the number of particles */
    nstep = atoi(argv[2]); /*command line input expected which provides the number of particles */


    /*initialize the position of particles in space and also their velocity and acceleration vector values*/
    initParticles(npart);
    clock_gettime(CLOCK_REALTIME, &bgn);
    /*begin the simulation with the initialization done previously and using simulation parameters defined in the global constants*/
    begin_simulation();

    clock_gettime(CLOCK_REALTIME, &nd);
    printf("%f\n", timediff(bgn,nd));
}
