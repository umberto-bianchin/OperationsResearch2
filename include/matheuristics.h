#ifndef MATHEURISTICS_H  
#define MATHEURISTICS_H

#include <cplex.h>
#include <tsp_utilities.h>
#include <utils.h>
#include <mincut.h>
#include <cplex_utilities.h>

static const double probabilities[4] = {0.2, 0.5, 0.8, 0.9};

void hard_fixing(instance *inst);
void set_warmup_solution(CPXENVptr env, CPXLPptr lp, instance *inst, solution *s);
void fix_random_edges(CPXENVptr env, CPXLPptr lp, instance *inst, double *xstar, double P);
void reset_lb(CPXENVptr env, CPXLPptr lp, instance *inst);

#endif /* MATHEURISTICS_H */