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
#define EPS_COST 	        1e-2  	// epsilon used to compare double costs
#define INF_COST 	        10e38  	// infinity for cost, used to represent infinity cost

#define PARAMS              17      // number of params that can be passed by command line
// -- List of indeces for the parameters array --
#define KICK                0       // kick in the instance params array
#define K_OPT               1       // value for the kopt: 3 for 3-opt, 5 for 5-opt, 6 or more for random k-opt
#define ALPHA               2       // alpha for grasp algorithm
#define MIN_COSTS           3       // number of minimum cost to store and choose randomly in grasp algorithm

#define MAX_TENURE          4       // max_tenure in the instance params
#define MIN_TENURE          5       // min_tenure in the instance params
#define TENURE_STEP         6       // tenure_step in the instance params
#define WARMUP              7       // use warmup solution with CPLEX
#define POSTING             8       // use warmup solution with CPLEX
#define CONCORDE            9       // use concorde with CPLEX
#define DEPTH               10      // posting solution for nodes <= depth
#define FIXEDPROB           11      // fixed probability to use for hard-fixing
#define PROBABILITY         12      // 1 = fixed probability, 0 decreasing for hard-fixing
#define K_LOCAL_BRANCHING   13      // local branching parameter
#define CDEPTH              14      // cplex fixing depth
#define POPULATION_SIZE     15      // population size for genetic algorithm
#define GENERATION_SIZE     16      // generation increment for genetic algorithm

#define ALGORITHMS_SIZE       10
static const char *algorithms[ALGORITHMS_SIZE] = {
    // If you add an algorithm here remember to add the corresponding choose_run_algorithm
    "N = Nearest Neighbour", 
    "E = Extra Mileage", 
    "V = Variable Neighbourhood Search", 
    "G = GRASP", 
    "T = Tabu Search",
    "B = Benders",
    "C = Branch and Cut",
    "H = Hard Fixing",
    "L = Local Branching",
    "P = Genetic Algorithm"
};

static const char *parameters[PARAMS] = {
    "kick = Kick parameter for VNS", 
    "kopt = K parameter for K-opt algorithm, [kopt >= 3]",

    "alpha = Alpha parameter for GRASP", 
    "minc = Min Costs parameter for GRASP", 
    
    "maxt = Max Tenure parameter for Tabu Search", 
    "mint = Min Tenure parameter for Tabu Search",
    "stept = Tenure Step parameter for Tabu Search",

    "warmup = Use warmup solution with CPLEX, [0] for NO, [1] for YES"
    "posting = Post a solution found by patching heuristic for Benders, [0] for NO, [1] for YES"
    "concorde = Use concorde for Branch and Cut, [0] for NO, [1] for YES"
    "depth = Posting solution only for nodes with depth <= depth"

    "fixedprob = [1] Hard fixing with fixed probability, [0] decreasing probability"
    "probability = Probability to use in hard fixing, integer probability (20 for 20 percent)"
    "klocal = k-opt neighborhood used in local branching"
    "cdepth = Depth until CPLEX runs for hard and soft fixing"

    "population = Population size for genetic algorithm"
    "generation = Generation increment for genetic algorithm"
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