#ifndef TSP_UTILITIES_H_  
#define TSP_UTILITIES_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <cplex.h>  
#include <stdbool.h>

#define VERBOSE		50
#define EPS_COST 	10e-5  // epsilon for cost, used to compare two double costs (instead of using ==)
#define INF_COST 	10e38  // infinity for cost, used to represent infinity cost

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
	int seed;
	double timelimit;
    char input_file[1000];

    //global data
	int *best_sol;
	double best_cost; 	
    
	int *solution;			// contains the current solution
	double solution_cost;

    double *cost;			// actual cost of the edges
    
	double time_limit;		// time limit in seconds
    double t_start;

    int integer_costs;

} instance;     

void print_error(const char *err);
void free_instance(instance *inst);
void choose_rand_sol(instance *inst);
void plot_solution(instance *inst);

void compute_all_costs(instance *inst);
bool check_solution(instance *inst);
void update_best_sol(instance *inst);
void refine_opt(instance *inst);

double dist(int i, int j, instance *inst);

#endif   /* TSP_UTILITIES_H_ */ 