#ifndef HEURISTICS_H_  
#define HEURISTICS_H_

#include <chrono.h>
#include <tsp_utilities.h>

void nearest_neighbour(instance *inst, int start_node);
void multi_start_nearest_neighbours(instance *inst, double timelimit);
void variable_neighbourhood(instance *inst, double timelimit);
void extra_mileage(instance *inst);
void grasp(instance *inst, int start_node);
void multi_start_grasp(instance *inst, double timelimit);
void tabu(instance *inst, double timelimit);
void genetic_algorithm(instance *inst, double timelimit);
void extra_mileage_operation(instance *inst, solution *s, int nInserted, int *inserted);
int short_cut(instance *inst, solution *s);
void generate_random_solutions(instance *inst, solution *solutions);
int generate_childs(instance *inst, solution *solutions, solution *childs, double timelimit, double best_cost);
void kill_members(instance *inst, solution *solutions, solution *childs, double total_cost, int best_idx, int generatedChild);

#endif   /* HEURISTICS_H_ */ 