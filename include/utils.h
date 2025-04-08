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
#define EPS_ERR 	        1e-5  	// epsilon used to compare if two double costs are equal (instead of using ==)
#define EPS_COST 	        300  	// epsilon used to compare double costs
#define INF_COST 	        10e38  	// infinity for cost, used to represent infinity cost
#define PARAMS              7       // number of params that can be passed by command line

// -- List of indeces for the parameters array --
#define KICK                0       // kick in the instance params array
#define K_OPT               1       // value for the kopt: 3 for 3-opt, 5 for 5-opt, 6 or more for random k-opt
#define ALPHA               2       // alpha for grasp algorithm
#define MIN_COSTS           3       // number of minimum cost to store and choose randomly in grasp algorithm

#define MAX_TENURE          4       // max_tenure in the instance params
#define MIN_TENURE          5       // min_tenure in the instance params
#define TENURE_STEP         6       // tenure_step in the instance params

#define ALGORITHMS_SIZE       6
static const char *algorithms[ALGORITHMS_SIZE] = {
    // If you add an algorithm here remember to add the corresponding choose_run_algorithm
    "N = Nearest Neighbour", 
    "E = Extra Mileage", 
    "V = Variable Neighbourhood Search", 
    "G = GRASP", 
    "T = Tabu Search",
    "B = Benders"
};

static const char *parameters[PARAMS] = {
    "kick = Kick parameter for VNS", 
    "kopt = K parameter for K-opt algorithm, [kopt >= 3]",

    "alpha = Alpha parameter for GRASP", 
    "minc = Min Costs parameter for GRASP", 
    
    "maxt = Max Tenure parameter for Tabu Search", 
    "mint = Min Tenure parameter for Tabu Search",
    "stept = Tenure Step parameter for Tabu Search"
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
void add_solution(solutions *sol, double cost, double time);         
void free_solution_struct(solutions *sol);              
void plot_solutions(instance *inst);
void plot_cplex_solutions(instance *inst);
void benchmark_algorithm_by_params(instance *inst);
void setAlgorithmId(instance *inst, char *algorithmID);

#endif /* UTILS_H_ */