#include <matheuristics.h>

void hard_fixing(instance *inst){
    printf("Initializing solution with heuristics\n");

	double initialization_timelimit = inst->time_limit/10.0;
	variable_neighbourhood(inst, initialization_timelimit);
	solution s;
	copy_solution(&s, &inst->best_solution, inst->nnodes);

    printf("Best solution found with heuristics: %f\n", s.cost);
	printf("Time taken: %f\n", second() - inst->t_start);

	double remaining_time = inst->time_limit - (second() - inst->t_start);
	double local_time_limit;
	int error = 0;
	int *indices = malloc(inst->nnodes * sizeof(int));
	char *lu = malloc(inst->nnodes * sizeof(char));
	double *bd = malloc(inst->nnodes * sizeof(double));

	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if (error) print_error("CPXcreateprob() error");

	build_model(inst ,env, lp);
	
	int cnt = 0;
	while(remaining_time > 0){
		// randon value between inst->params[MIN_HARD] and inst->params[MAX_HARD]				
		double P = inst->params[MIN_HARD] + ((double)rand() / (double)RAND_MAX )* (inst->params[MAX_HARD] - inst->params[MIN_HARD]);
		P /= 100.0; // convert to percentage

		for(int i = 0; i < inst->nnodes + 1; i++){
			if(((double)rand() / (double)RAND_MAX) < P){
				indices[cnt] = s.path[i];
				lu[cnt] = 'L'; 	// lower bound
				bd[cnt++] = 1.0; 
			}
		}

		// https://www.ibm.com/docs/en/cofz/12.9.0?topic=cpxxchgbds-cpxchgbds
		error = CPXchgbds(env, lp, cnt, indices, lu, bd);
		if (error){
			printf("CPX error code %d\n", error);
			print_error("hard_fixing() error"); 
		}

		// error = CPXmipopt(env, lp);
		// Update solution s with the new solution from CPLEX

		
		if(VERBOSE >= INFO){
			printf("Iter %4d, lower bound %10.2lf, ncomp %4d, time %5.2lf\n", cnt, s.cost, cnt, second() - inst->t_start);
		}

		cnt = 0;	
		local_time_limit = remaining_time/10.0;
		remaining_time = inst->time_limit - (second() - inst->t_start);
	}

	free_solution(&s);
}
