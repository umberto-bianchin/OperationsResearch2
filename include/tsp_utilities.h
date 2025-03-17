#ifndef TSP_UTILITIES_H_  
#define TSP_UTILITIES_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <cplex.h>
#include <chrono.h>
#include <stdbool.h>
#include <ctype.h>

/**
 * @brief
 * Struct used to store the best solutions found during the search dinamically
 */
typedef struct solution_struct{   
    double *all_best_cost;
    int size;
    int capacity;
} solutions;


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

	solutions history_best_cost;	// contains all the best costs found during the search
} instance;



void initialize_instance(instance *inst);
void allocate_instance(instance *inst);
void copy_instance(instance *inst, instance *new_inst);
void free_instance(instance *inst);

void choose_rand_sol(instance *inst);
void check_solution(instance *inst, bool best);
void update_best_solution(instance *inst);
double compute_solution_cost(instance *inst, int *tour);

void compute_all_costs(instance *inst);
double calculate_delta(int i, int j, instance *inst);
void reverse_segment(int start, int end, instance *inst);
double dist(int i, int j, instance *inst);

void two_opt(instance *inst);
int find_best_move(instance *inst, int a, int b, int c, int d, int e, int f, int n);
void apply_best_move(instance *inst, int i, int j, int k, int best_case);
void three_opt(instance *inst);

#endif   /* TSP_UTILITIES_H_ */ 