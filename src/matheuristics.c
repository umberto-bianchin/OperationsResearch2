#include <matheuristics.h>

void hard_fixing(instance *inst){
	int error = 0;
	CPXENVptr env = CPXopenCPLEX(&error); // Initialize CPLEX environment
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); // Create CPLEX problem instance
	if (error) print_error("CPXcreateprob() error");

	CPXLONG contextid = CPX_CALLBACKCONTEXT_RELAXATION | CPX_CALLBACKCONTEXT_CANDIDATE;

	build_model(inst, env, lp); // Build the optimization model
	
	// CPLEX's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF); // Disable screen output
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-9); // Set optimality gap
	if(VERBOSE >= DEBUG) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Enable screen output if verbose

	if(CPXcallbacksetfunc(env, lp, contextid, cpx_callback, inst)){
		print_error("CPXcallbacksetfunc() error");
	}
	
	int *comp = (int *) calloc(inst->nnodes, sizeof(int)); // Component array for connected components
	int *succ = (int *) calloc(inst->nnodes, sizeof(int)); // Successor array for solution
	int ncomp = 9999;
	inst->ncols = CPXgetnumcols(env, lp); // Get number of variables
	
	printf("Initializing solution with heuristics\n");

	double initialization_timelimit = inst->time_limit/10.0;
	nearest_neighbour(inst, rand() % inst->nnodes); // Initialize solution using nearest neighbor heuristic
	solution s;
	copy_solution(&s, &inst->best_solution, inst->nnodes);
	two_opt(inst, &s, inst->time_limit); // Improve solution using 2-opt heuristic
	variable_neighbourhood(inst, initialization_timelimit); // Apply variable neighborhood search
	
	copy_solution(&s, &inst->best_solution, inst->nnodes);

	int max_edges = inst->ncols;
	int *index = (int *) calloc(max_edges, sizeof(int)); // Array for variable indices
	double *xstar = (double *) calloc(max_edges, sizeof(double)); // Array for variable values

	solution_to_CPX(&s, inst->nnodes, index, xstar); // Convert solution to CPLEX format

	printf("Best solution found with heuristics: %f, after time %f\n", s.cost, second() - inst->t_start);

	double remaining_time = inst->time_limit - (second() - inst->t_start);
	double local_time_limit;
	int iteration = 0;
	double new_cost = CPX_INFBOUND;
	double old_cost = s.cost;
	srand(inst->seed); // Seed random number generator

	while(remaining_time > 0){
		set_warmup_solution(env, lp, inst, &s); // Set warm-up solution for CPLEX

		double P;
		if(iteration < 4){
			P = probabilities[3]; // High probability of fixing edges
		} else if (iteration < 6){
			P = probabilities[2];
		} else if (iteration < 8) {
			P = probabilities[1];
		} else {
			P = probabilities[0]; // Low probability of fixing edges
		}

		fix_random_edges(env, lp, inst, xstar, P); // Fix random edges based on probability

		local_time_limit = (inst->time_limit/10.0 < remaining_time) ? inst->time_limit/10.0 : remaining_time;
		CPXsetdblparam(env, CPX_PARAM_TILIM, local_time_limit); // Set time limit for optimization

		error = CPXmipopt(env,lp); // Solve the MIP problem

		if (error){
			print_error("CPXmipopt() error"); 
		}
		
		error = CPXgetobjval(env, lp, &new_cost); // Get objective value of the solution
		if (error) {
			print_error("CPXgetobjval() error");
		}

		if(new_cost < old_cost){ // Update solution if a better one is found
			old_cost = new_cost;
			if (CPXgetx(env, lp, xstar, 0, inst->ncols-1)){
				print_error("CPXgetx() error");
			}
	
			build_sol(xstar, inst, succ, comp, &ncomp); // Build solution from CPLEX variables

			solution_from_CPX(inst, &s, succ); // Update solution structure
		}

		reset_lb(env, lp, inst); // Reset lower bounds of variables

		printf("Best solution found at iteration %3d: %6.4f, after time %f\n", iteration, old_cost, second() - inst->t_start);
		iteration++;
		remaining_time = inst->time_limit - (second() - inst->t_start);
	}

	update_best_solution(inst, &s); // Update the best solution in the instance
	free(xstar);
	free(index);
	free_route(&s);
	free(succ);
	free(comp);
	CPXfreeprob(env, &lp); // Free CPLEX problem
	CPXcloseCPLEX(&env); // Close CPLEX environment
}

void set_warmup_solution(CPXENVptr env, CPXLPptr lp, instance *inst, solution *s){
	int *index = (int *) calloc(inst->ncols, sizeof(int)); // Array for variable indices
	double *xstar = (double *) calloc(inst->ncols, sizeof(double)); // Array for variable values

	solution_to_CPX(s, inst->nnodes, index, xstar); // Convert solution to CPLEX format

	int effortlevel = CPX_MIPSTART_NOCHECK; // Effort level for warm-up solution
	int beg = 0;
	int error = CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, index, xstar, &effortlevel, NULL); // Add warm-up solution
	if (error) {
		print_error("CPXaddmipstarts() error");
	}

	free(xstar);
	free(index);
}

void fix_random_edges(CPXENVptr env, CPXLPptr lp, instance *inst, double *xstar, double P){
	for(int i = 0; i < inst->nnodes - 1; i++){
		for (int j = i + 1; j < inst->nnodes; j++){
			int varidx = xpos(i, j, inst->nnodes); // Get variable index for edge (i, j)
			if (i != j && xstar[varidx] > 0.5 && ((double)rand() / RAND_MAX) < P){ // Fix edge with probability P
				char lu = 'L';
				double bd = 1.0;

				int error = CPXchgbds(env, lp, 1, &varidx, &lu, &bd); // Change lower bound to fix edge
				if(error){
					print_error("CPXchgbds() error");
				}
				
			}
		}
	}
}

void reset_lb(CPXENVptr env, CPXLPptr lp, instance *inst){
	for(int i = 0; i < inst->nnodes - 1; i++){
		for (int j = i + 1; j < inst->nnodes; j++){
			int varidx = xpos(i, j, inst->nnodes); // Get variable index for edge (i, j)
			char lu = 'L';
			double bd = 0.0;

			int error = CPXchgbds(env, lp, 1, &varidx, &lu, &bd); // Reset lower bound to 0
			if(error){
					print_error("CPXchgbds() error");
			}
		}			
	}
}