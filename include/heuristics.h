#ifndef HEURISTICS_H_  
#define HEURISTICS_H_

#include <tsp_utilities.h>
#include <chrono.h>

void nearest_neighbour(instance *inst, int start_node);
void all_nearest_neighbours(instance *inst);
void vns(instance *inst);
void extra_mileage(instance *inst);
void grasp(instance *inst, double alpha);

#endif   /* HEURISTICS_H_ */ 