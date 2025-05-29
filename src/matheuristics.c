#include <matheuristics.h>

/**
 * @brief
 * CPLEX Fixing Matheuristic:
 * - Hard fixing (algorithm == 'H'):
 *     Iteratively fixes a random subset of incumbent edges (lower‐bound = 1)
 *     with decreasing probability, then reoptimizes to improve the solution.
 * - Local branching (algorithm == 'L'):
 *     Adds a cut to restrict the search to the neighborhood of the incumbent,
 *     by requiring that at least (n − k) of the incumbent’s edges remain in the tour.
 *
 * This function builds the TSP model, generates an initial heuristic solution
 * (Nearest Neighbor + 2-opt + VNS), then loops until time_limit:
 *   1. Injects the current best solution as a warm start.
 *   2. Applies either hard fixing or local branching based on inst->algorithm.
 *   3. Solves the MIP for a fixed time slice.
 *   4. Updates the incumbent if an improvement is found.
 *
 * @param inst
 *   Pointer to a fully initialized TSP instance, with time_limit and params set.
 */
void cplex_fixing(instance *inst){
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
	
	CPXsetintparam(env, CPX_PARAM_NODELIM, inst->params[CDEPTH]);
	//CPXsetintparam(env, CPX_PARAM_MIPEMPHASIS, CPX_MIPEMPHASIS_FEASIBILITY);	

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
	if(inst->algorithm == 'H' && inst->params[FIXEDPROB]){
		if(inst->params[PROBABILITY] > MAX_PROB){
			printf("Probability %d too high, fixing it to %d percent\n", inst->params[PROBABILITY], MAX_PROB);
			P = MAX_PROB / 100;
		} else {
			P = inst->params[PROBABILITY] / 100;
		}
	}
	
	if(inst->algorithm == 'L'){
		inst->params[CDEPTH] = 9000;	// best tuning for local branching algorithm
	}
	
	while(remaining_time > 0){
		// Probability not fixed: start high, then decrease
		if(inst->algorithm == 'H' && !inst->params[FIXEDPROB]){
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

		if(inst->algorithm == 'H'){
			fix_random_edges(env, lp, inst, xstar, P);
		} else {
			add_local_branching(env, lp, inst, xstar);
		}

		// Allocate time slice for MIP
		//local_time_limit = (inst->time_limit/10.0 < remaining_time) ? inst->time_limit/10.0 : remaining_time;
		CPXsetdblparam(env, CPX_PARAM_TILIM, remaining_time);
		
		error = CPXmipopt(env,lp);

		if (error){
			print_error("CPXmipopt() error"); 
		}
		
		error = CPXgetobjval(env, lp, &new_cost);
		if (error) {
			print_error("CPXgetobjval() error");
		}
		
		add_solution(&(inst->history_best_costs), inst->best_solution.cost, -1);

		if(new_cost < old_cost - EPS_COST){
			old_cost = new_cost;
			if (CPXgetx(env, lp, xstar, 0, inst->ncols-1)){
				print_error("CPXgetx() error");
			}
	
			build_sol(xstar, inst, succ, comp, &ncomp);

			solution_from_CPX(inst, &s, succ);
		}

		if(inst->algorithm == 'H'){
			// Reset lower bounds of variables
			reset_lb(env, lp, inst); 
		} else {
			// Remove the last constraint added
			int last_row = CPXgetnumrows(env, lp) - 1;
			CPXdelrows(env, lp, last_row, last_row);
		}

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

void add_local_branching(CPXENVptr env, CPXLPptr lp, instance *inst, double *xstar){
	// Local branching cut
	// n - k elements have to be fixed
	int k = inst->params[K_LOCAL_BRANCHING];
	int n = inst->nnodes;
	int n_fixed = n - k;

	int rmatbeg = 0;
	char sense = 'G';
	double rhs = (double) n_fixed;
	int indices[n];
	double coefs[n];

	for (int i = 0; i < n; ++i){
		int a = inst->best_solution.path[i];
		int b = inst->best_solution.path[(i+1)%n];
		int h = xpos(a, b, inst->nnodes);

		indices[i] = h;
		coefs[i] = 1.0;
	}
	if (CPXaddrows(env, lp, 0, 1, n, &rhs, &sense, &rmatbeg, indices, coefs, NULL, NULL)){
		print_error("CPXaddrows() error");
	}
}