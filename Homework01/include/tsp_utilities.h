#ifndef TSP_UTILITIES_H_  
#define TSP_UTILITIES_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <cplex.h>  

#define VERBOSE		50

/**
 * @brief 
 * TSP data structure
 */
typedef struct {   
	
	//input data
	int nnodes;
	double *xcoord;
	double *ycoord;

	// parameters 
	int randomseed;
	double timelimit;
    char input_file[1000];

    //global data
	double *best_sol;

} instance;     


void print_error(const char *err);
void free_instance(instance *inst);
void choose_rand_sol(instance *inst);
void plot_solution(instance *inst);

#endif   /* TSP_UTILITIES_H_ */ 