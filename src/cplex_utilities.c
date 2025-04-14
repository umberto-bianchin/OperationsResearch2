#define _CRT_SECURE_NO_WARNINGS
#include <cplex_utilities.h>

/**
 * @brief 
 * Function that calculates the position of edge [i, j] in the CPLEX TSP model
 * @param i the first node
 * @param j the second node
 * @param inst the tsp instance
 * @return int the position of the edge in the CPLEX model
 */
int xpos(int i, int j, instance *inst){
	if (i == j) print_error(" i == j in xpos" );
	if (i > j) return xpos(j,i,inst);
	int pos = i * inst->nnodes + j - ((i + 1) * (i + 2)) / 2;
	return pos;
}

/**
 * @brief 
 * Build the CPLEX model adding all the variables and all the constraints needed
 * @param inst the tsp instance
 * @param env the CPLEX environment variable
 * @param lp the CPLEX problem
 */
void build_model(instance *inst, CPXENVptr env, CPXLPptr lp){    
	int izero = 0; 
	char binary = 'B'; 
	char **cname = (char **) calloc(1, sizeof(char *));
	cname[0] = (char *) calloc(100, sizeof(char));

	// Add binary variables x(i,j) for i < j  
	for(int i = 0; i < inst->nnodes; i++){
		for(int j = i+1; j < inst->nnodes; j++){
			sprintf(cname[0], "x(%d,%d)", i+1,j+1);
			double obj = dist(i, j, inst);
			double lb = 0.0;
			double ub = 1.0;
            
			if (CPXnewcols(env, lp, 1, &obj, &lb, &ub, &binary, cname)) print_error(" wrong CPXnewcols on x var.s");
    		if (CPXgetnumcols(env,lp)-1 != xpos(i, j, inst)) print_error(" wrong position for x var.s");
		}
	} 

    // Add the degree constraints 
	int *index = (int *) calloc(inst->nnodes, sizeof(int));
	double *value = (double *) calloc(inst->nnodes, sizeof(double));
	double rhs = 2.0;
	char sense = 'E';

	for (int h = 0; h < inst->nnodes; h++){
		sprintf(cname[0], "degree(%d)", h+1);
        int nnz = 0;

		for (int i = 0; i < inst->nnodes; i++){
			if (i == h) continue;
			
			index[nnz] = xpos(i, h, inst);
			value[nnz] = 1.0;
			nnz++;
		}

		if (CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]))
			print_error("CPXaddrows(): error 1");
	} 

	free(value);
	free(index);

	free(cname[0]);
	free(cname);  
}

/**
 * @brief 
 * Call the function to build CPLEX model and compute the best solution using CPLEX
 * @param inst the tsp instance
 * @return int 0 if no error
 */
int TSPopt(instance *inst){  
	// Open CPLEX model
	int error = 0;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if (error) print_error("CPXcreateprob() error");

	build_model(inst, env, lp);
	
	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	//CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->seed);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-9);
	if(VERBOSE >= DEBUG) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen

	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int ncols = CPXgetnumcols(env, lp);
	double *xstar = (double *) calloc(ncols, sizeof(double));
	int ncomp = 9999;
	double objval = 0.0;
	int iter = 0;

	while(ncomp >= 2){
		double elapsed_time = second() - inst->t_start;
		double residual_time = inst->time_limit - elapsed_time;
		CPXsetdblparam(env, CPX_PARAM_TILIM, residual_time);

		error = CPXmipopt(env,lp);

		elapsed_time = second() - inst->t_start;
		residual_time = inst->time_limit - elapsed_time;

		if(residual_time <= 0){
			printf("Exceded time limit while computing CPXmipopt(), exiting the loop\n");
			break;
		}
	
		if (error){
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}

		error = CPXgetobjval(env, lp, &objval);
		if (error) print_error("CPXgetbestobjval() error");
	
		add_solution(&(inst->history_best_costs), objval, elapsed_time);

		if (CPXgetx(env, lp, xstar, 0, ncols-1))
			print_error("CPXgetx() error");
		
		build_sol(xstar, inst, succ, comp, &ncomp);

		if(VERBOSE >= INFO){
			printf("Iter %4d, lower bound %10.2lf, ncomp %4d, time %5.2lf\n", iter, objval, ncomp, second() - inst->t_start);
			fflush(NULL);
		}

		if(ncomp >= 2){
			add_sec(inst, env, lp, comp, ncomp, ncols);
		}
		
		iter++;
	}	

	// Write the model in an appropriate file 
	if (VERBOSE >= DEBUG)
		CPXwriteprob(env, lp, "history/model.lp", NULL);

	
	if(ncomp >= 2){
		if(VERBOSE >= INFO){
			printf("Entering Patching Heuristic method with ncomp %4d, time %5.2lf\n", ncomp, second() - inst->t_start);
			fflush(NULL);
		}

		patching_heuristic(inst, succ, comp, &ncomp);
		objval = -1.0;

		if(VERBOSE >= INFO){
			printf("Exiting Patching Heuristic method with ncomp %4d, time %5.2lf\n", ncomp, second() - inst->t_start);
			fflush(NULL);
		}

	}

	copy_best_solution(inst, env, lp, succ, objval);

	// Free and close cplex model   
	free(xstar);
	free(succ);
	free(comp);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	return 0;
}

// Build succ() and comp() wrt xstar()...
void build_sol(const double *xstar, instance *inst, int *succ, int *comp, int *ncomp){

	// Needed only for debug purposes
	#ifdef DEBUGON
		int *degree = (int *) calloc(inst->nnodes, sizeof(int));
		
		for (int i = 0; i < inst->nnodes; i++){
			for (int j = i+1; j < inst->nnodes; j++){
				int k = xpos(i,j,inst);

				if (fabs(xstar[k]) > EPS_ERR && fabs(xstar[k]-1.0) > EPS_ERR )
					print_error(" wrong xstar in build_sol()");
				
				if (xstar[k] > 0.5){
					++degree[i];
					++degree[j];
				}
			}
		}

		for ( int i = 0; i < inst->nnodes; i++ ){
			if (degree[i] != 2)
				print_error("wrong degree in build_sol()");
		}	
		free(degree);
	#endif

	*ncomp = 0;
	for (int i = 0; i < inst->nnodes; i++){
		succ[i] = -1;
		comp[i] = -1;
	}
	
	for (int start = 0; start < inst->nnodes; start++){
		
		// Node "start" was already visited, skip it
		if (comp[start] >= 0) continue;

		// New component is found
		(*ncomp)++;
		int i = start;
		int done = 0;

		// Visit the current component
		while (!done){
			comp[i] = *ncomp;
			done = 1;

			for (int j = 0; j < inst->nnodes; j++){

				// The edge [i,j] is selected in xstar and j was not visited before
				if (i != j && xstar[xpos(i,j,inst)] > 0.5 && comp[j] == -1){
					succ[i] = j;
					i = j;
					done = 0;
					break;
				}
			}
		}	
		succ[i] = start;
	}
}

void add_sec(instance *inst, CPXENVptr env, CPXLPptr lp, int *comp, int ncomp, int ncols){
	int izero = 0; 
	char **cname = (char **) calloc(1, sizeof(char *));
	int max_edges = ncols * (ncols - 1) / 2;
	int *index = (int *) calloc(max_edges, sizeof(int));
	double *value = (double *) calloc(max_edges, sizeof(double));
	cname[0] = (char *) calloc(100, sizeof(char));

	for(int k = 1; k <= ncomp; k++){
		double rhs = -1;
		char sense = 'L';
		int nnz = 0;

		for(int i=0; i < inst->nnodes; i++){
			if(comp[i] != k)
				continue;

			rhs++;		
			for(int j=i+1; j < inst->nnodes; j++){
				if(comp[j] != k) continue;

				value[nnz] = 1.0;
				index[nnz] = xpos(i, j, inst);
				nnz++;
			}
		}

		int size = rhs + 1.1;
		if(nnz != size * (size-1)/2) {
			printf("nnz %d, size %d, expected nnz %d\n", nnz, size, (size * (size-1)/2));
			print_error("Wrong nnz in add_sec()");
		}

		sprintf(cname[0], "sec(%d)", k);
		CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]);
	}

	free(value);
	free(index);
	free(cname[0]);
	free(cname);
}

/**
 * @brief 
 * Copy the optimal solution found by CPLEX into inst->best_solution
 */
void copy_best_solution(instance *inst, CPXENVptr env, CPXLPptr lp, int *succ, double objval) {
	solution s;
    if (inst->best_solution.path == NULL) {
        inst->best_solution.path = (int *) calloc(inst->nnodes, sizeof(int));
    }

	s.path = (int *) calloc(inst->nnodes, sizeof(int));
    
    int current = 0; 
    for (int i = 0; i < inst->nnodes + 1; i++) {
		s.path[i] = current;
        current = succ[current];
    }

	double total_cost = 0.0;
	for (int i = 0; i < inst->nnodes; i++)
		total_cost += inst->costs[s.path[i] * inst->nnodes + s.path[i +1]];
	
	if(objval != -1.0){
		if(fabs(objval - total_cost) > EPS_ERR){
			free_route(&s);
			free(inst);
			print_error("Calculated cost is different from CPLEX cost");
		}
	}

	s.cost = total_cost;
	check_solution(inst, &s);

	inst->best_solution.path = s.path;
	inst->best_solution.cost = total_cost;
	add_solution(&(inst->history_best_costs), total_cost, second() - inst->t_start);
}

/**
 * @brief 
 * Method that merges different component given by CPLEX algorithm
 * @param inst the tsp instance
 * @param succ the array of successors founded by CPLEX that will be modified
 * @param comp the array of components founded by CPLEX that will be modified
 * @param ncomp the number of component that will be modified
 */
void patching_heuristic(instance *inst, int *succ, int *comp, int *ncomp){
    if (*ncomp < 2){
        return;
    }

	int merged_components = 0;
	for (int c1 = 1; c1 <= *ncomp; c1++){
		for (int c2 = c1 + 1; c2 <= *ncomp; c2++){
			int *nodes_c1 = (int *) calloc(inst->nnodes, sizeof(int *));
			int *nodes_c2 = (int *) calloc(inst->nnodes, sizeof(int *));
			int size_c1 = 0, size_c2 = 0;
			int best_i1 = -1, best_i2 = -1, best_j1 = -1, best_j2 = -1;
			double min_delta = INF_COST;
			bool orientation;

			for (int i = 0; i < inst->nnodes; i++){
				if (comp[i] == c1){
					nodes_c1[size_c1++] = i;
				}
				else if(comp[i] == c2) {
					nodes_c2[size_c2++] = i;
				}
			}

			for (int i = 0; i < size_c1; i++) {
				for (int j = 0; j < size_c2; j++) {
					int i1 = nodes_c1[i], j1 = nodes_c2[j];
					int i2 = succ[i1], j2 = succ[j1];

					double delta1 = delta_cost(inst, i1, j1, i2, j2, true);
					double delta2 = delta_cost(inst, i1, j1, i2, j2, false);

					if (delta1 < delta2 && delta1 < min_delta) {
						min_delta = delta1;
						best_i1 = i1;
						best_i2 = i2;
						best_j1 = j1;
						best_j2 = j2;
						orientation = true;
					} else if (delta2 < min_delta){
						min_delta = delta2;
						best_i1 = i1;
						best_i2 = i2;
						best_j1 = j1;
						best_j2 = j2;
						orientation = false;
					}
				}
			}

			if (best_i1 != -1) {
				if(orientation){
					reverse_cycle(inst, best_j1, succ);

					succ[best_i1] = best_j1;
					succ[best_j2] = best_i2;
				}
				else{
					succ[best_i1] = best_j2;
					succ[best_j1] = best_i2;
				}

				for(int j = 0; j < size_c2; j++)
					comp[nodes_c2[j]] = c1;

				merged_components++;
			}
		}
	}

	*ncomp -= merged_components;
}

/**
 * @brief 
 * Calculate the cost of creating two new edges when two components are merged
 * @param inst the tsp instance
 * @param i1 the first node in the first component
 * @param j1 the first node in the second component
 * @param i2 the second node in the first component
 * @param j2 the second node in the second component
 * @param option boolean to decide how to merge the two components
 * @return double 
 */
double delta_cost(instance *inst, int i1, int j1, int i2, int j2, bool option){
	if(option)
    	return (inst->costs[i1 * inst->nnodes + j1] + inst->costs[i2 * inst->nnodes + j2]) - (inst->costs[i1 * inst->nnodes + i2] + inst->costs[j1 * inst->nnodes + j2]);
	else
		return (inst->costs[i1 * inst->nnodes + j2] + inst->costs[i2 * inst->nnodes + j1]) - (inst->costs[i1 * inst->nnodes + i2] + inst->costs[j1 * inst->nnodes + j2]);	
}

/**
 * @brief 
 * Function to reversing a component cycle
 * @param inst the tsp instance
 * @param start the starting node
 * @param succ the successors array founded by CPLEX
 */
void reverse_cycle(instance *inst, int start, int *succ){
    int current = start;
    int prev = -1;
    
    do {
        int next = succ[current];
        succ[current] = prev;
        prev = current;
        current = next;
    } while (current != start);

	succ[start] = prev;
}