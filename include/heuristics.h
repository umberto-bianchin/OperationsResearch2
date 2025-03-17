#ifndef HEURISTICS_H_  
#define HEURISTICS_H_

#include <chrono.h>
#include <tsp_utilities.h>

void nearest_neighbour(instance *inst, int start_node);
void multi_start_nearest_neighbours(instance *inst);
void variable_neighbourhood(instance *inst);
void extra_mileage(instance *inst);
void grasp(instance *inst, int start_node);
void multi_start_grasp(instance *inst);

#endif   /* HEURISTICS_H_ */ 