#ifndef MATHEURISTICS_H  
#define MATHEURISTICS_H

#include <cplex.h>
#include <tsp_utilities.h>
#include <utils.h>
#include <mincut.h>
#include <cplex_utilities.h>

#define MAX_PROB 90

static const double probabilities[4] = {0.2, 0.4, 0.5, 0.8};

void cplex_fixing(instance *inst);
void fix_random_edges(CPXENVptr env, CPXLPptr lp, instance *inst, double *xstar, double P);
void reset_lb(CPXENVptr env, CPXLPptr lp, instance *inst);
void add_local_branching(CPXENVptr env, CPXLPptr lp, instance *inst, double *xstar);

#endif /* MATHEURISTICS_H */