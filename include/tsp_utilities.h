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
typedef struct solutions_struct{   
    double *all_costs;
    int size;
    int capacity;
} solutions;


typedef struct solution_struct{   
	int *path;
	double cost;
} solution;

/**
 * @brief 
 * TSP data structure
 */
typedef struct instance_struct{   
	// input data
	int nnodes;
	double *xcoord;
	double *ycoord;
	char algorithm;
	char running_mode;
	int *params;

	// parameters 
	int seed;
    char input_file[1000];

	solution best_solution;		// best route found

    double *costs;			// actual costs of the edges
    
	double time_limit;		// time limit in seconds
    double t_start;

	solutions history_best_costs;	// contains all the best costs found during the search
	solutions history_costs;			// contains all the costs found during the search
} instance;

void initialize_instance(instance *inst);
void allocate_instance(instance *inst);
void allocate_route(solution *s, int nnodes);
void copy_instance(instance *inst, instance *new_inst);
void free_route(solution *s);
void free_instance(instance *inst);
void set_random_coord(instance *inst);
void copy_solution(solution *dest, solution *src, int nnodes);
void choose_rand_sol(instance *inst);
void check_solution(instance *inst, solution *s);
void update_best_solution(instance *inst, solution *s);
void compute_solution_cost(instance *inst, solution *s);

void compute_all_costs(instance *inst);
double calculate_delta(int i, int j, instance *inst, solution *s);
void reverse_segment(int start, int end, solution *s);
double dist(int i, int j, instance *inst);

void two_opt(instance *inst, solution *s);
int find_best_move(instance *inst, int a, int b, int c, int d, int e, int f, int n);
void apply_best_move(instance *inst, int i, int j, int k, int best_case, solution *s);
void three_opt(instance *inst, solution *s);

#endif   /* TSP_UTILITIES_H_ */ 