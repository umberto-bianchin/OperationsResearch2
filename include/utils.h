#ifndef UTILS_H_  
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <tsp_utilities.h>
#include <heuristics.h>

#define INFO		        20
#define ERROR		        50
#define DEBUG		        100
#define VERBOSE		        INFO
#define MAX_NO_IMPROVEMENT  10000   //number of iterations in which variable neihbourhood can get no improvement
#define KICK                5       //number of time that 3-opt is called in variable neighbourhood
#define EPS_COST 	        10e-5  	// epsilon for cost, used to compare two double costs (instead of using ==)
#define INF_COST 	        10e38  	// infinity for cost, used to represent infinity cost


void print_error(const char *err, bool terminate);
void plot_solution(instance *inst, bool best);
void choose_run_algorithm(instance *inst);
void benchmark_algorithm_by_time(instance *inst);

#endif /* UTILS_H_ */ 