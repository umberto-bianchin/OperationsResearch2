#ifndef MATHEURISTICS_H  
#define MATHEURISTICS_H

#include <cplex.h>
#include <tsp_utilities.h>
#include <utils.h>
#include <mincut.h>
#include <cplex_utilities.h>

static const double probabilities[4] = {0.2, 0.5, 0.8, 0.9};

void hard_fixing(instance *inst);

#endif /* MATHEURISTICS_H */