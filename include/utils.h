#ifndef UTILS_H_  
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <tsp_utilities.h>
#include <heuristics.h>
#include <csv_parser.h>

#define INFO		        20
#define ERROR		        50
#define DEBUG		        100
#define VERBOSE		        INFO
#define MAX_COORDINATES     10000    // used in random nodes generator, it's the upperbound of the x and y values

#define MAX_NO_IMPROVEMENT  5000    // number of iterations in which variable neihbourhood can get no improvement
#define KICK                5       // number of time that 3-opt is called in variable neighbourhood
#define EPS_COST 	        10e-5  	// epsilon for cost, used to compare two double costs (instead of using ==)
#define INF_COST 	        10e38  	// infinity for cost, used to represent infinity cost

#define ALPHA               0.2     // alpha for grasp algorithm
#define MIN_COSTS             3     // number of minimum cost to store and choose randomly in grasp algorithm

#define MAX_TENURE            500    // size of the tabu list in tabu search
#define MIN_TENURE            100    // size of the tabu list in tabu search
#define TENURE_STEP           50     // size of the tabu list in tabu search

#define ALGORITHMS_SIZE       5
static const char *algorithms[ALGORITHMS_SIZE] = {
    // If you add an algorithm here remember to add the corresponding choose_run_algorithm
    "N = Nearest Neighbour", 
    "E = Extra Mileage", 
    "V = Variable Neighbourhood Search", 
    "G = GRASP", 
    "T = Tabu Search"
};

void check_valid_algorithm(char algorithm);
const char* print_algorithm(char algorithm);
void print_algorithms();

void print_error(const char *err);
void plot_solution(instance *inst, solution *s);
void choose_run_algorithm(instance *inst);
void benchmark_algorithm_by_time(instance *inst);
void allocate_solution_struct(solutions *sol);          
void add_solution(solutions *sol, double cost);         
void free_solution_struct(solutions *sol);              
void plot_solutions(instance *inst);
void benchmark_algorithm_by_params(instance *inst);

#endif /* UTILS_H_ */