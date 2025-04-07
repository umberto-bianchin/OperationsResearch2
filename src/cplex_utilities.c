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
	CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->seed);	
	if(VERBOSE >= DEBUG) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen

	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int ncomp = 9999;
	double *xstar;
	double objval = 0;

	while(ncomp >= 2){
		double elapsed_time = (second() - inst->t_start);
		double residual_time = inst->time_limit - elapsed_time;
		CPXsetdblparam(env, CPX_PARAM_TILIM, residual_time);

		error = CPXmipopt(env,lp);
	
		if (error){
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}

		CPXgetbestobjval(env, lp, &objval);
		if (error) print_error("CPXgetbestobjval() error");
	
		add_solution(&(inst->history_best_costs), objval, elapsed_time);

		// Use the optimal solution found by CPLEX
		int ncols = CPXgetnumcols(env, lp);
		xstar = (double *) calloc(ncols, sizeof(double));

		if (CPXgetx(env, lp, xstar, 0, ncols-1))
			print_error("CPXgetx() error");
		
		build_sol(xstar, inst, succ, comp, &ncomp);

		add_sec(inst, env, lp, comp, &ncomp);
	}	

	inst->best_solution.cost = objval;

	/*for (int i = 0; i < inst->nnodes; i++){
		for (int j = i+1; j < inst->nnodes; j++){
			if ( xstar[xpos(i,j,inst)] > 0.5 ) printf("  ... x(%3d,%3d) = 1\n", i+1,j+1);
		}
	}*/

	// Write the model in an appropriate file 
	if (VERBOSE >= DEBUG)
		CPXwriteprob(env, lp, "history/model.lp", NULL); 

	// Free and close cplex model   
	free(xstar);
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

void add_sec(instance *inst, CPXENVptr env, CPXLPptr lp, int *comp, int *ncomp){
	int izero = 0; 
	char **cname = (char **) calloc(1, sizeof(char *));

	for(int k = 1; k <= (*ncomp); k++){
		int max_edges = inst->nnodes * (inst->nnodes - 1);
		int *index = (int *) calloc(max_edges, sizeof(int));
		double *value = (double *) calloc(max_edges, sizeof(double));
		cname[0] = (char *) calloc(100, sizeof(char));

		double rhs = -1;
		char sense = 'L';
		int nnz = 0;

		for(int i=0; i < inst->nnodes; i++){
			if(comp[i] != k)
				continue;

			rhs++;		
			for(int j=0; j < inst->nnodes; j++){
				if(j==i) continue;
				if(comp[j] != k) continue;

				value[nnz] = 1.0;
				index[nnz] = xpos(i, j, inst);
				nnz++;
			}
		}
		
		sprintf(cname[0], "sec(%d)", k);
		CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0]);
		
		(*ncomp)--;
		
		free(index);
		free(value);
		free(cname[0]);
	}

	free(cname);
}