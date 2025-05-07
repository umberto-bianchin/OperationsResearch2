#include <matheuristics.h>

void hard_fixing(instance *inst){
	int error = 0;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if (error) print_error("CPXcreateprob() error");

	CPXLONG contextid = CPX_CALLBACKCONTEXT_RELAXATION | CPX_CALLBACKCONTEXT_CANDIDATE;

	build_model(inst, env, lp);
	
	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-9);
	if(VERBOSE >= DEBUG) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen

	if(CPXcallbacksetfunc(env, lp, contextid, cpx_callback, inst)){
		print_error("CPXcallbacksetfunc() error");
	}
	
	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int ncomp = 9999;
	inst->ncols = CPXgetnumcols(env, lp);
	
    printf("Initializing solution with heuristics\n");

	double initialization_timelimit = inst->time_limit/10.0;
	nearest_neighbour(inst, rand() % inst->nnodes);     //need to initialize and optimize the first solution
	solution s;
	copy_solution(&s, &inst->best_solution, inst->nnodes);
	two_opt(inst, &s, inst->time_limit);
	variable_neighbourhood(inst, initialization_timelimit);
	
	copy_solution(&s, &inst->best_solution, inst->nnodes);

	int max_edges = inst->ncols;
	int *index = (int *) calloc(max_edges, sizeof(int));
	double *xstar = (double *) calloc(max_edges, sizeof(double));

	printf("CIAO nodi %d, costo %f\n", inst->nnodes, s.cost);
	solution_to_CPX(&s, inst->nnodes, index, xstar);
	printf("CIAO 2\n");

	/*int effortlevel = CPX_MIPSTART_NOCHECK;
	int beg = 0;
    // Set the warm start solution in CPLEX
    int error = CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, index, xstar, &effortlevel, NULL);
    if (error) {
        print_error("CPXaddmipstarts() error");
    }*/

    printf("Best solution found with heuristics: %f, after time %f\n", s.cost, second() - inst->t_start);

	double remaining_time = inst->time_limit - (second() - inst->t_start);
	double local_time_limit;
	int iteration = 0;
	double new_cost = CPX_INFBOUND;
	double old_cost = s.cost;
	srand(inst->seed);

	while(remaining_time > 0){		
		double P;
		if(iteration < 4){
			P = probabilities[3];
		} else if (iteration < 6){
			P = probabilities[2];
		} else if (iteration < 8) {
			P = probabilities[1];
		} else {
			P = probabilities[0];
		}
		// 0.2 e 0.9 li metterei come parametri in ingresso
		// minimo random e massimo random 				
		//double P = 0.2 + ((double)rand() / (double)RAND_MAX )* (0.9 - 0.2);

		for(int i = 0; i < inst->nnodes - 1; i++){
			for (int j = i + 1; j < inst->nnodes; j++){
				int varidx = xpos(i, j, inst->nnodes); 
				if (i != j && xstar[varidx] > 0.5 && ((double)rand() / RAND_MAX) < P){
					char lu = 'L';                       // 3) sense of bound
            		double bd = 1.0;

					error = CPXchgbds(env, lp, 1, &varidx, &lu, &bd);
					if(error){
						print_error("CPXchgbds() error");
					}
					
				}
			}
		}

		local_time_limit = (inst->time_limit/20.0 < remaining_time) ? inst->time_limit/20.0 : remaining_time;
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

		for(int i = 0; i < inst->nnodes - 1; i++){
			for (int j = i + 1; j < inst->nnodes; j++){
				int varidx = xpos(i, j, inst->nnodes);
				char lu = 'L';
				double bd = 0.0;

				CPXchgbds(env, lp, 1, &varidx, &lu, &bd);
			}			
		}

		printf("Best solution found at iteration %3d: %6.4f, after time %f\n", iteration, old_cost, second() - inst->t_start);
		iteration++;
		remaining_time = inst->time_limit - (second() - inst->t_start);
	}

	free(xstar);
	free(index);
	free_route(&s);
	free(succ);
	free(comp);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env);
}
