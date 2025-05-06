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

	while(remaining_time > 0){
		// 0.2 e 0.9 li metterei come parametri in ingresso
		// minimo random e massimo random 				
		double P = 0.2 + ((double)rand() / (double)RAND_MAX )* (0.9 - 0.2);

		for(int i = 0; i < inst->nnodes + 1; i++){
			if(((double)rand() / (double)RAND_MAX) < P){
				
			}
		}

		local_time_limit = remaining_time/10.0;

		remaining_time = inst->time_limit - (second() - inst->t_start);
	}

	free_solution(&s);
}
