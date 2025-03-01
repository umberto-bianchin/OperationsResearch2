#ifndef TSP_H_  

#define TSP_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 

#include <cplex.h>  

#define VERBOSE				    50		// printing level  (=10 only incumbent, =20 little output, =50-60 good, =70 verbose, >=100 cplex log)                             

//data structures  

typedef struct {   
	
	//input data
	int nnodes;
	double *xcoord;
	double *ycoord;

	// parameters 
	int randomseed;
	double timelimit;						// overall time limit, in sec.s
    char input_file[1000];		  			// input file

    //global data
	double *best_sol;						// best sol. available

} instance;        

#endif   /* TSP_H_ */ 
