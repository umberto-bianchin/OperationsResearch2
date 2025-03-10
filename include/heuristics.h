#ifndef HEURISTICS_H_  
#define HEURISTICS_H_

#include <tsp_utilities.h>
#include <chrono.h>

void nearest_neighbour(instance *inst, int start_node);
void all_nearest_neighbours(instance *inst);
void extra_mileage(instance *inst);

#endif   /* HEURISTICS_H_ */ 