#ifndef TSP_UTILITIES_H_  
#define TSP_UTILITIES_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <cplex.h>
#include <chrono.h>
#include <utils.h>

#define EPS_COST 	10e-5  	// epsilon for cost, used to compare two double costs (instead of using ==)
#define INF_COST 	10e38  	// infinity for cost, used to represent infinity cost

typedef struct{
	int *tour;
	double cost;
} solution;

/**
 * @brief 
 * TSP data structure
 */
typedef struct {   
	
	// input data
	int nnodes;
	double *xcoord;
	double *ycoord;

	// parameters 
	int seed;
    char input_file[1000];

	solution best_solution;

    double *costs;			// actual costs of the edges
    
	double time_limit;		// time limit in seconds
    double t_start;

    int integer_costs;

} instance;     

void initialize_instance(instance *inst);
void allocate_instance(instance *inst);

void free_instance(instance *inst);
void choose_rand_sol(instance *inst);
void plot_solution(instance *inst, char best);
void choose_run_algorithm(instance *inst);

void compute_all_costs(instance *inst);
void check_solution(instance *inst, char best);
void update_best_solution(instance *inst);
void calc_solution_cost(instance *inst);
double calculate_delta(int i, int j, instance *inst);
void swap_nodes(int i, int j, instance *inst);
void two_opt(instance *inst);

double dist(int i, int j, instance *inst);

#endif   /* TSP_UTILITIES_H_ */ 