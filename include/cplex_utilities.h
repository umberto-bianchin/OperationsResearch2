#ifndef CPLEX_UTILITIES_H_  
#define CPLEX_UTILITIES_H_

#include <cplex.h>
#include <tsp_utilities.h>
#include <utils.h>

#define DEBUGON    // Comment out to avoid debugging 

int xpos(int i, int j, instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
int TSPopt(instance *inst);
void add_sec(instance *inst, CPXENVptr env, CPXLPptr lp, int *comp, int *ncomp);
void copy_best_solution(instance *inst, CPXENVptr env, CPXLPptr lp, int *succ, int *comp, int *ncomp);

#endif   /*CPLEX_UTILITIES_H_ */ 