#ifndef UTILS_H_  
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include <tsp_utilities.h>
#include <heuristics.h>
#include <csv_parser.h>

#define INFO		        20
#define ERROR		        50
#define DEBUG		        100
#define VERBOSE		        INFO
#define MAX_COORDINATES     10000    // used in random nodes generator, it's the upperbound of the x and y values

#define MAX_NO_IMPROVEMENT  5000    // number of iterations in which variable neihbourhood can get no improvement
#define EPS_ERR 	        10e-5  	// epsilon used to compare if two double costs are equal (instead of using ==)
#define EPS_COST 	        100  	// epsilon used to compare double costs
#define INF_COST 	        10e38  	// infinity for cost, used to represent infinity cost
#define PARAMS              6       // number of params that can be passed by command line

#define KICK                0       // index of the kick in the instance params array
#define ALPHA               1       // index of the alpha for grasp algorithm
#define MIN_COSTS           2       // index of the number of minimum cost to store and choose randomly in grasp algorithm

#define MAX_TENURE          3       // index of the max_tenure in the instance params
#define MIN_TENURE          4       // index of the min_tenure in the instance params
#define TENURE_STEP         5       // index of the tenure_step in the instance params

#define ALGORITHMS_SIZE       5
static const char *algorithms[ALGORITHMS_SIZE] = {
    // If you add an algorithm here remember to add the corresponding choose_run_algorithm
    "N = Nearest Neighbour", 
    "E = Extra Mileage", 
    "V = Variable Neighbourhood Search", 
    "G = GRASP", 
    "T = Tabu Search"
};

static const char *parameters[PARAMS] = {
    "kick = Kick param for VNS", 
    "alpha = Alpha param for GRASP", 
    "minc = Min Costs param for GRASP", 
    "maxt = Max Tenure param for Tabu Search", 
    "mint = Min Tenure param for Tabu Search",
    "stept = Tenure Step param for Tabu Search"
};

void check_valid_algorithm(char algorithm);
const char* print_algorithm(char algorithm);
void print_algorithms();
void print_parameters();

void print_error(const char *err);
void plot_solution(instance *inst, solution *s);
void choose_run_algorithm(instance *inst);
void benchmark_algorithm_by_time(instance *inst);
void allocate_solution_struct(solutions *sol);          
void add_solution(solutions *sol, double cost);         
void free_solution_struct(solutions *sol);              
void plot_solutions(instance *inst);
void benchmark_algorithm_by_params(instance *inst);
void setAlgorithmId(instance *inst, char *algorithmID);

#endif /* UTILS_H_ */