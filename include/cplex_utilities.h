#ifndef CPLEX_UTILITIES_H_  
#define CPLEX_UTILITIES_H_

#include <cplex.h>
#include <tsp_utilities.h>
#include <utils.h>
#include <mincut.h>

//#define DEBUGON    // Comment out to avoid debugging 

typedef struct callback_params{
	CPXCALLBACKCONTEXTptr context;
	instance* inst;
} relaxation_callback_params;

static int CPXPUBLIC cpx_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);
static int CPXPUBLIC sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);
static int CPXPUBLIC concorde_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle);
static int violated_cuts_callback(double cutval, int n, int *members, void *userhandle);

int xpos(int i, int j, instance *inst);
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp);
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp);
int TSPopt(instance *inst);
void add_sec(instance *inst, CPXENVptr env, CPXLPptr lp, CPXCALLBACKCONTEXTptr context, CPXLONG contextid, int *comp, int ncomp, bool callback);
void copy_best_solution(instance *inst, solution *s, int *succ, double objval);
void patching_heuristic(instance *inst, int *succ, int *comp, int *ncomp, double *objval);
double delta_cost(instance *inst, int i1, int j1, int i2, int j2, bool option);
void reverse_cycle(instance *inst, int start, int *succ);
void warmup_CPX_solution(instance *inst, CPXENVptr env, CPXLPptr lp);
void solution_to_CPX(instance *inst, int *index, double *value);
void post_CPX_solution(instance *inst, CPXCALLBACKCONTEXTptr context, int *succ, int *comp, int *ncomp, double *objval);


#endif   /*CPLEX_UTILITIES_H_ */ 