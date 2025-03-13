#ifndef TSP_UTILITIES_H_  
#define TSP_UTILITIES_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <cplex.h>
#include <chrono.h>
#include <utils.h>
#include <ctype.h>

#define EPS_COST 	10e-5  	// epsilon for cost, used to compare two double costs (instead of using ==)
#define INF_COST 	10e38  	// infinity for cost, used to represent infinity cost

/**
 * @brief 
 * TSP data structure
 */
typedef struct instance_struct{   
	
	// input data
	int nnodes;
	double *xcoord;
	double *ycoord;

	// parameters 
	int seed;
    char input_file[1000];

	int *best_solution;
	double best_cost;

	int *solution;
	double solution_cost;

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
void check_solution(instance *inst, bool best);
void update_best_solution(instance *inst);
double compute_solution_cost(instance *inst, int *tour);
double calculate_delta(int i, int j, instance *inst);
void reverse_segment(int start, int end, instance *inst);
void two_opt(instance *inst);
double find_best_move(instance *inst, int a, int b, int c, int d, int e, int f, int n, double *totCost);
void apply_best_move(instance *inst, int i, int j, int k, int best_case);
void three_opt(instance *inst);
double dist(int i, int j, instance *inst);

#endif   /* TSP_UTILITIES_H_ */ 