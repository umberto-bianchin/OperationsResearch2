#include <matheuristics.h>

/**
 * @brief
 * Hard-fixing matheuristic: iteratively fixes a subset of edges and solves to improve solution.
 *
 * Builds model, initializes heuristic solution (NN + 2-opt + VNS), then loops:
 * - injects warm start
 * - randomly fixes edges with decreasing probability
 * - solves MIP for a time slice
 * - updates best solution
 *
 * @param inst  Pointer to initialized TSP instance with time_limit and params set.
 */
void hard_fixing(instance *inst){
	int error = 0;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1");
	if (error) print_error("CPXcreateprob() error");

	CPXLONG contextid = CPX_CALLBACKCONTEXT_RELAXATION | CPX_CALLBACKCONTEXT_CANDIDATE;

	build_model(inst, env, lp);
	
	// CPLEX's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-9);
	if(VERBOSE >= DEBUG) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON);

	// Register callbacks for SEC separation
	if(CPXcallbacksetfunc(env, lp, contextid, cpx_callback, inst)){
		print_error("CPXcallbacksetfunc() error");
	}
	
	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int ncomp = 9999;
	inst->ncols = CPXgetnumcols(env, lp);

	warmup_CPX_solution(inst, env, lp, true);
	
	int max_edges = inst->ncols;
	int *index = (int *) calloc(max_edges, sizeof(int));
	double *xstar = (double *) calloc(max_edges, sizeof(double));

	solution_to_CPX(&inst->best_solution, inst->nnodes, index, xstar);

	// Main hard-fixing loop
	double remaining_time = inst->time_limit - (second() - inst->t_start);
	double local_time_limit;
	int iteration = 0;
	double new_cost = CPX_INFBOUND;
	double old_cost = inst->best_solution.cost;
	double P;
	solution s;
	allocate_route(&s, inst->nnodes);
	copy_solution(&s, &inst->best_solution, inst->nnodes);	// setting first solution as the solution found with heuristic
	srand(inst->seed); // ensure repeatability

	// Setting fixed probability, with maximum value MAX_PROB
	if(inst->params[FIXEDPROB]){
		if(inst->params[PROBABILITY] > MAX_PROB){
			printf("Probability %d too high, fixing it to %d percent\n", inst->params[PROBABILITY], MAX_PROB);
			P = MAX_PROB / 100;
		} else {
			P = inst->params[PROBABILITY] / 100;
		}
	}
	
	while(remaining_time > 0){
		// Probability not fixed: start high, then decrease
		if(!inst->params[FIXEDPROB]){
			if(iteration < 4){
			P = probabilities[3];
			} else if (iteration < 6){
				P = probabilities[2];
			} else if (iteration < 8) {
				P = probabilities[1];
			} else {
				P = probabilities[0];
			}
		} 

		fix_random_edges(env, lp, inst, xstar, P);

		// Allocate time slice for MIP
		local_time_limit = (inst->time_limit/10.0 < remaining_time) ? inst->time_limit/10.0 : remaining_time;
		CPXsetdblparam(env, CPX_PARAM_TILIM, local_time_limit);

		error = CPXmipopt(env,lp);

		if (error){
			print_error("CPXmipopt() error"); 
		}
		
		error = CPXgetobjval(env, lp, &new_cost);
		if (error) {
			print_error("CPXgetobjval() error");
		}

		if(new_cost < old_cost){
			old_cost = new_cost;
			if (CPXgetx(env, lp, xstar, 0, inst->ncols-1)){
				print_error("CPXgetx() error");
			}
	
			build_sol(xstar, inst, succ, comp, &ncomp);

			solution_from_CPX(inst, &s, succ);
		}

		// Reset lower bounds of variables
		reset_lb(env, lp, inst); 

		printf("Best solution found at iteration %3d: %6.4f, after time %f\n", iteration, old_cost, second() - inst->t_start);
		iteration++;
		remaining_time = inst->time_limit - (second() - inst->t_start);
		set_warmup_solution(env, lp, inst, &s); // Set warm-up solution for CPLEX for the next iteration
	}

	update_best_solution(inst, &s); 
	free(xstar);
	free(index);
	free_route(&s);
	free(succ);
	free(comp);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
}

/**
 * @brief
 * Randomly fix a subset of tour edges by setting lower bounds to 1.
 *
 * @param env   CPLEX environment pointer.
 * @param lp    CPLEX problem pointer.
 * @param inst  TSP instance with nnodes.
 * @param xstar Binary solution vector (1 if edge in tour).
 * @param P     Probability to fix each present edge.
 */
void fix_random_edges(CPXENVptr env, CPXLPptr lp, instance *inst, double *xstar, double P){
	for(int i = 0; i < inst->nnodes - 1; i++){
		for (int j = i + 1; j < inst->nnodes; j++){
			int varidx = xpos(i, j, inst->nnodes); 
			if (i != j && xstar[varidx] > 0.5 && ((double)rand() / RAND_MAX) < P){ 
				char lu = 'L';
				double bd = 1.0;

				int error = CPXchgbds(env, lp, 1, &varidx, &lu, &bd);
				if(error){
					print_error("CPXchgbds() error");
				}
				
			}
		}
	}
}

/**
 * @brief
 * Reset all edge lower bounds back to 0 for next iteration.
 *
 * @param env  CPLEX environment pointer.
 * @param lp   CPLEX problem pointer.
 * @param inst TSP instance with nnodes.
 */
void reset_lb(CPXENVptr env, CPXLPptr lp, instance *inst){
	for(int i = 0; i < inst->nnodes - 1; i++){
		for (int j = i + 1; j < inst->nnodes; j++){
			int varidx = xpos(i, j, inst->nnodes);
			char lu = 'L';
			double bd = 0.0;

			int error = CPXchgbds(env, lp, 1, &varidx, &lu, &bd);
			if(error){
					print_error("CPXchgbds() error");
			}
		}			
	}
}