#ifndef HEURISTICS_H_  
#define HEURISTICS_H_

#include <tsp_utilities.h>
#include <chrono.h>

#define KICK    5

void nearest_neighbour(instance *inst, int start_node);
void multi_start_nearest_neighbours(instance *inst);
void variable_neighbourhood(instance *inst);
void extra_mileage(instance *inst);
void grasp(instance *inst, double alpha);

#endif   /* HEURISTICS_H_ */ 