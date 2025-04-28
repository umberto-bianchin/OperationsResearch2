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
int xpos(int i, int j, int nnodes){
	if (i == j) print_error(" i == j in xpos" );
	if (i > j) return xpos(j,i,nnodes);
	int pos = i * nnodes + j - ((i + 1) * (i + 2)) / 2;

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
    		if (CPXgetnumcols(env,lp)-1 != xpos(i, j, inst->nnodes)) print_error(" wrong position for x var.s");
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
			
			index[nnz] = xpos(i, h, inst->nnodes);
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
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-9);
	if(VERBOSE >= DEBUG) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	
	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	inst->ncols = CPXgetnumcols(env, lp);
	double *xstar = (double *) calloc(inst->ncols, sizeof(double));
	int ncomp = 9999;
	double objval = 0.0;
	int iter = 0;
	
	if(inst->params[WARMUP]){
		warmup_CPX_solution(inst, env, lp);
	}

	if(inst->algorithm == 'C'){
		CPXsetintparam(env, CPX_PARAM_THREADS, 1); 	// just for debugging
		CPXLONG contextid;

		if(inst->params[CONCORDE]){
			contextid = CPX_CALLBACKCONTEXT_RELAXATION | CPX_CALLBACKCONTEXT_CANDIDATE;
		} else{
			contextid = CPX_CALLBACKCONTEXT_CANDIDATE;
		}
		if(CPXcallbacksetfunc(env, lp, contextid, cpx_callback, inst)){
			print_error("CPXcallbacksetfunc() error");
		}
	}

	while(ncomp >= 2){
		double elapsed_time = second() - inst->t_start;
		double residual_time = inst->time_limit - elapsed_time;
		CPXsetdblparam(env, CPX_PARAM_TILIM, residual_time);
		
		error = CPXmipopt(env,lp);

		if (error){
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
		}

		elapsed_time = second() - inst->t_start;
		residual_time = inst->time_limit - elapsed_time;

		if(residual_time <= 0){
			printf("Exceded time limit while computing CPXmipopt(), exiting the loop\n");
			break;
		}
		
		error = CPXgetobjval(env, lp, &objval);
		if (error) {
			char errmsg[CPXMESSAGEBUFSIZE];
			CPXgeterrorstring(NULL, error, errmsg);
			print_error(errmsg);
		}

		if(inst->algorithm == 'B'){
			add_solution(&(inst->history_best_costs), objval, elapsed_time);
		}
		
		if (CPXgetx(env, lp, xstar, 0, inst->ncols-1)){
			print_error("CPXgetx() error");
		}

		build_sol(xstar, inst, succ, comp, &ncomp);

		if(VERBOSE >= INFO){
			if(inst->algorithm == 'B'){
				printf("Iter %4d, lower bound %10.2lf, ncomp %4d, time %5.2lf\n", iter, objval, ncomp, second() - inst->t_start);
			}
			else{
				printf("Final obj %10.2lf, ncomp %4d, time %5.2lf\n", objval, ncomp, second() - inst->t_start);
			}
			fflush(NULL);
		}

		if(inst->algorithm == 'B' && ncomp >= 2){
			add_sec(inst, env, lp, NULL, -1, comp, ncomp, false);
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

		patching_heuristic(inst, succ, comp, &ncomp, &objval);

		if(VERBOSE >= INFO){
			printf("Exiting Patching Heuristic method with ncomp %4d, time %5.2lf\n", ncomp, second() - inst->t_start);
			fflush(NULL);
		}
	}

	solution s;
	allocate_route(&s, inst->nnodes);
	copy_best_solution(inst, &s, succ, objval);

	// Free and close cplex model
	free_route(&s);
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
				int k = xpos(i,j,inst->nnodes);

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
				if (i != j && xstar[xpos(i,j,inst->nnodes)] > 0.5 && comp[j] == -1){
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

/**
 * @brief 
 * Driver method used to decide whether we are in Branch and Cut or Concorde Branch and Cut
 * @param context CPLEX callback context
 * @param contextid The callback context identifier
 * @param userhandle Pointer to the user-defined data structure (the TSP instance)
 * @return 0 if successful, 1 otherwise
 */
static int CPXPUBLIC cpx_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle){ 
	instance* inst = (instance*) userhandle;
	if(contextid == CPX_CALLBACKCONTEXT_CANDIDATE) return sec_callback(context, contextid, inst); 
	if(contextid == CPX_CALLBACKCONTEXT_RELAXATION) return concorde_callback(context, contextid, inst);
	print_error("Contextid unknownn in cpx_callback()"); 
	return 1; 	
}

/**
 * @brief 
 * Method used in the Branch and Cut execution without concorde
 * @param context CPLEX callback context
 * @param contextid The callback context identifier (must be CPX_CALLBACKCONTEXT_CANDIDATE)
 * @param userhandle Pointer to the user-defined data structure (the TSP instance)
 * @return 0 if successful, 1 otherwise
 */
static int CPXPUBLIC sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle){ 
	instance* inst = (instance*) userhandle;  
	double* xstar = (double*) malloc(inst->ncols * sizeof(double));  
	double objval = CPX_INFBOUND;

	if (CPXcallbackgetcandidatepoint(context, xstar, 0, inst->ncols-1, &objval)){
		print_error("CPXcallbackgetcandidatepoint error");
	}

	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	int ncomp = 9999;

	build_sol(xstar, inst, succ, comp, &ncomp);

	if(ncomp >= 2){
		add_sec(inst, NULL, NULL, context, contextid, comp, ncomp, true);
	}
	
	double elapsed_time = second() - inst->t_start;
	add_solution(&inst->history_best_costs, objval, elapsed_time);

	if(inst->params[POSTING]){
		CPXLONG nodedepth;
		int status;
	
		int error = CPXcallbackgetinfoint(context, CPXCALLBACKINFO_CANDIDATE_SOURCE, &status);

		if(error){
			print_error("CPXgetcallbackinfo() error getting CANDIDATE_SOURCE");
		}

		if(status == CPX_LAZYCONSTRAINTCALLBACK_NODE){
			error = CPXcallbackgetinfolong(context, CPXCALLBACKINFO_NODEDEPTH, &nodedepth);
			if(error) {
				print_error("CPXgetcallbackinfo() error getting NODEDEPTH");
			}
		}
		
		if (nodedepth <= inst->params[DEPTH]){
			post_CPX_solution(inst, context, succ, comp, &ncomp, &objval);
		}
	}
	
	free(succ);
	free(comp);
	free(xstar);
	return 0; 
}

/**
 * @brief 
 * Method used in the Concorde Branch and Cut execution
 * @param context CPLEX callback context
 * @param contextid The callback context identifier (must be CPX_CALLBACKCONTEXT_RELAXATION)
 * @param userhandle Pointer to the user-defined data structure (the TSP instance)
 * @return 0 if successful, 1 otherwise
 */
static int CPXPUBLIC concorde_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle){
	instance* inst = (instance*) userhandle;  
	double* xstar = (double*) malloc(inst->ncols * sizeof(double));  
	double objval = CPX_INFBOUND;

	if (CPXcallbackgetrelaxationpoint(context, xstar, 0, inst->ncols-1, &objval)){
		print_error("CPXcallbackgetrelaxationpoint error");
	}

	int* elist =(int *) calloc(2 * inst->ncols, sizeof(int)); // elist contains each pair of vertex such as (1,2), (1,3), (1,4), (2, 3), (2,4), (3,4) so in list becomes: 1,2,1,3,1,4,2,3,2,4,3,4
	int *comps = NULL;
    int *compscount = NULL;
    int num_edges = 0, k = 0, ncomp = 0;

	for (int i = 0; i < inst->nnodes; i++){
		for (int j = i + 1; j < inst->nnodes; j++){
			elist[k++] = i;
			elist[k++] = j;
			num_edges++;
		}
	}

	// Checking whether or not the graph is connected with the fractional solution.
	if (CCcut_connect_components(inst->nnodes, num_edges, elist, xstar, &ncomp, &compscount, &comps)){
		print_error("CCcut_connect_components() error ");
	}

	if (ncomp == 1){
		relaxation_callback_params params = { .context = context, .inst = inst };

		if (CCcut_violated_cuts(inst->nnodes, inst->ncols, elist, xstar, 2.0 - EPS_ERR, violated_cuts_callback, &params)){
			print_error("%d, CCcut_violated_cuts() error ");
		}
	} else {
		int *comp = (int *) calloc(inst->nnodes, sizeof(int));
		int *succ = (int *) calloc(inst->nnodes, sizeof(int));
		int ncomp = 9999;

		build_sol(xstar, inst, succ, comp, &ncomp);
	
		add_sec(inst, NULL, NULL, context, contextid, comp, ncomp, true);

		free(succ);
		free(comp);
	}


    free(compscount);
    free(comps);
	free(elist);
    free(xstar);
	return 0;
}

/**
 * @brief
 * Callback function for handling fractional solutions in the CPLEX Branch-and-Cut framework,
 * leveraging Concorde's mincut routines to detect and separate violated SEC constraints.
 * @param context CPLEX callback context
 * @param contextid The callback context identifier (must be CPX_CALLBACKCONTEXT_RELAXATION)
 * @param userhandle Pointer to the user-defined data structure (the TSP instance)
 * @return 0 if successful, 1 otherwise
 */
static int violated_cuts_callback(double cutval, int n, int *members, void *userhandle){
    relaxation_callback_params *params = (relaxation_callback_params *)userhandle;
    instance *inst = params->inst;
    CPXCALLBACKCONTEXTptr context = params->context;

    double rhs = n - 1;
    char sense = 'L';
    int matbeg = 0;

    int *index = (int *) malloc(sizeof(int) * (n * (n - 1)) / 2);
    double *value = (double *) malloc(sizeof(double) * (n * (n - 1)) / 2);

    int nnz = 0;
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            int u = members[i];
            int v = members[j];
            index[nnz] = xpos(u, v, inst->nnodes);
            value[nnz] = 1.0;
            nnz++;
        }
    }

    int purgeable = CPX_USECUT_FILTER;
    int local = 0;
    int status = CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &matbeg, index, value, &purgeable, &local);

    free(index);
    free(value);
    return status;
}

/**
 * @brief
 * Adds or rejects Subtour Elimination Constraints (SECs) based on component connectivity.
 * If `callback` is true, the SECs are added via:
 * - `CPXcallbackrejectcandidate()` in the candidate context
 * - `CPXcallbackaddusercuts()` in the relaxation context
 * If `callback` is false, the SECs are added directly to the model using `CPXaddrows()`.
 *
 * @param inst The TSP instance
 * @param env The CPLEX environment (used only when not in a callback)
 * @param lp The CPLEX problem (used only when not in a callback)
 * @param context The callback context if inside a callback (NULL otherwise)
 * @param contextid The context ID (RELAXATION or CANDIDATE) used only if callback is true
 * @param comp An array indicating the component ID of each node
 * @param ncomp The number of connected components in the current solution
 * @param callback Boolean flag indicating whether we're in a callback context
 */
void add_sec(instance *inst, CPXENVptr env, CPXLPptr lp, CPXCALLBACKCONTEXTptr context, CPXLONG contextid, int *comp, int ncomp, bool callback){
	int izero = 0; 
	char **cname = (char **) calloc(1, sizeof(char *));
	int max_edges = inst->nnodes * (inst->nnodes - 1) / 2;

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
				index[nnz] = xpos(i, j, inst->nnodes);
				nnz++;
			}
		}

		int size = rhs + 1.1;
		if(nnz != size * (size-1)/2) {
			printf("nnz %d, size %d, expected nnz %d\n", nnz, size, (size * (size-1)/2));
			print_error("Wrong nnz in add_sec()");
		}

		int rows = CPXgetnumrows(env, lp);
		sprintf(cname[0], "sec(%d, %d)", k, rows);

		if(callback){
			if(contextid == CPX_CALLBACKCONTEXT_CANDIDATE){
				if(CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value)){
					print_error("CPXcallbackrejectcandidate() error");
				}
			} else if(contextid == CPX_CALLBACKCONTEXT_RELAXATION){
				int purgeable = CPX_USECUT_FILTER;
				int local = 0;

				if(CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &izero, index, value, &purgeable, &local)){
					print_error("CPXcallbackaddusercuts() error");
				}
			}
		} else{
			if(CPXaddrows(env, lp, 0, 1, nnz, &rhs, &sense, &izero, index, value, NULL, &cname[0])){
				print_error("CPXaddrows() error");
			}
		}
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
void copy_best_solution(instance *inst, solution *s, int *succ, double objval) {
    if (inst->best_solution.path == NULL) {
        inst->best_solution.path = (int *) calloc(inst->nnodes, sizeof(int));
    }

    int current = 0; 
    for (int i = 0; i <= inst->nnodes; i++) {
		s->path[i] = current;
        current = succ[current];
    }

	double total_cost = 0.0;
	for (int i = 0; i < inst->nnodes; i++)
		total_cost += inst->costs[s->path[i] * inst->nnodes + s->path[i +1]];
	
	if(objval != -1.0){
		if(fabs(objval - total_cost) > EPS_ERR){
			free_route(s);
			free(inst);
			print_error("Calculated cost is different from CPLEX cost");
		}
	}

	s->cost = total_cost;
	update_best_solution(inst, s);
	add_solution(&(inst->history_best_costs), total_cost, second() - inst->t_start);
}

/**
 * @brief 
 * Method that merges different components given by CPLEX algorithm, the final solution will be different from the initial one
 * @param inst the tsp instance
 * @param succ the array of successors founded by CPLEX that will be modified
 * @param comp the array of components founded by CPLEX that will be modified
 * @param ncomp the number of component that will be modified
 */
void patching_heuristic(instance *inst, int *succ, int *comp, int *ncomp, double *objval){
    if (*ncomp < 2){
        return;
    }

	int merged_components = 0;
	for (int c1 = 1; c1 <= *ncomp; c1++){
		for (int c2 = c1 + 1; c2 <= *ncomp; c2++){
			bool valid_c2 = false;
            for (int k = 0; k < inst->nnodes; k++) {
                if (comp[k] == c2) {
                    valid_c2 = true;
                    break;
                }
            }
            if (!valid_c2) continue;

			int *nodes_c1 = (int *) calloc(inst->nnodes, sizeof(int));
			int *nodes_c2 = (int *) calloc(inst->nnodes, sizeof(int));
			int size_c1 = 0, size_c2 = 0;
			int best_i = -1, succ_i = -1, best_j = -1, succ_j = -1;
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
						best_i = i1;
						succ_i = i2;
						best_j = j1;
						succ_j = j2;
						orientation = true;
					} else if (delta2 < min_delta){
						min_delta = delta2;
						best_i = i1;
						succ_i = i2;
						best_j = j1;
						succ_j = j2;
						orientation = false;
					}
				}
			}

			if(best_i == -1 || best_j == -1)
				print_error("No valid pair best_i1 and best_j1 founded");
			
			if(orientation){
				reverse_cycle(inst, best_j, succ);

				succ[best_i] = best_j;
				succ[succ_j] = succ_i;
			}
			else{
				succ[best_i] = succ_j;
				succ[best_j] = succ_i;
			}

			for(int j = 0; j < size_c2; j++)
				comp[nodes_c2[j]] = c1;

			merged_components++;
			
			free(nodes_c2);
			free(nodes_c1);
		}
	}

	*ncomp -= merged_components;
	*objval = -1.0;
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

/**
 * @brief
 * Sets up a warmup solution for CPLEX
 * @param inst the tsp instance containing the solution
 * @param env the CPLEX environment
 * @param lp the CPLEX problem
 * @return int 0 if no error occurred
 */
void warmup_CPX_solution(instance *inst, CPXENVptr env, CPXLPptr lp) {
	double initialization_timelimit = inst->time_limit/10.0;
	nearest_neighbour(inst, rand() % inst->nnodes);     //need to initialize and optimize the first solution
	solution s;
	copy_solution(&s, &inst->best_solution, inst->nnodes);
	two_opt(inst, &s, initialization_timelimit);

	double remaining_time = initialization_timelimit - (second() - inst->t_start);
	variable_neighbourhood(inst, remaining_time);

    // Check if we have a valid solution to use as warm start
    if (inst->best_solution.path == NULL) {
        print_error("No solution available in inst->best_solution.path\n");
    }

	int max_edges = inst->ncols;// * (inst->ncols - 1) / 2;

	int *index = (int *) calloc(max_edges, sizeof(int));
	double *value = (double *) calloc(max_edges, sizeof(double));

	solution_to_CPX(&s, inst->nnodes, index, value);
    
	int effortlevel = CPX_MIPSTART_NOCHECK;
	int beg = 0;
    // Set the warm start solution in CPLEX
    int error = CPXaddmipstarts(env, lp, 1, inst->ncols, &beg, index, value, &effortlevel, NULL);
    if (error) {
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(NULL, error, errmsg);
		fprintf(stderr, "CPLEX Error %d: %s\n", error, errmsg);
        print_error("Failed to add MIP start, CPXaddmipstarts failed");
    }
    printf("Initial solution found with cost %lf after %5.2lf seconds", inst->best_solution.cost, (second() - inst->t_start));
    free(value);
    free(index);
	free_route(&s);
}

/**
 * @brief 
 * Posts a feasible heuristic solution to CPLEX from a callback context.
 * @param inst The TSP instance
 * @param context The CPLEX callback context
 * @param succ The array of successor nodes in the current solution
 * @param comp The array of component identifiers
 * @param ncomp Pointer to the number of components
 * @param objval Pointer to the current objective value
 */
void post_CPX_solution(instance *inst, CPXCALLBACKCONTEXTptr context, int *succ, int *comp, int *ncomp, double *objval){
	double elapsed_time = second() - inst->t_start;
	double residual_time = (inst->time_limit) - elapsed_time;

	patching_heuristic(inst, succ, comp, ncomp, objval);
	solution s;
	allocate_route(&s, inst->nnodes);
	copy_best_solution(inst, &s, succ, *objval);

	two_opt(inst, &s, residual_time);


	int max_edges = inst->ncols;
	int *index = (int *) calloc(max_edges, sizeof(int));
	double *xstar = (double *) calloc(max_edges, sizeof(double));

	solution_to_CPX(&s, inst->nnodes, index, xstar);

	int error = 0;
	if(s.cost < inst->best_solution.cost - EPS_COST){
		error = CPXcallbackpostheursoln(context, max_edges, index, xstar, s.cost, CPXCALLBACKSOLUTION_NOCHECK);
		if(VERBOSE >= INFO && !error) printf("Posted solution of cost %lf (my incubement is %lf)\n", s.cost, inst->best_solution.cost);
	}
	if(error){
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(NULL, error, errmsg);
		fprintf(stderr, "CPLEX Error %d: %s\n", error, errmsg);
		print_error("CPXcallbackpostheursoln() error");
	}

	free(xstar);
	free(index);
	free_route(&s);
} 

/**
 * @brief Converts the TSP solution from path representation to CPLEX variable representation.
 * 
 * This function translates the solution stored in `inst->best_solution.path` into the format
 * required by CPLEX. It populates the `index` array with the indices of the variables
 * corresponding to the edges in the solution and the `xstar` array with the values of these
 * variables (1.0 for edges included in the solution, 0.0 otherwise).
 * 
 * @param inst The TSP instance containing the solution and problem data.
 * @param index Array to store the indices of the variables corresponding to the edges in the solution.
 * @param xstar Array to store the values of the variables (binary: 1.0 for included edges, 0.0 otherwise).
 */
void solution_to_CPX(solution *s, int nnodes, int *index, double *xstar){
    // Populate the `index` array with the indices of all possible edges in the model.
    int k = 0;
    for (int i = 0; i < nnodes - 1; i++) {
        for (int j = i + 1; j < nnodes; j++) {
            index[k++] = xpos(i, j, nnodes);
        }
    }

    // Populate the `xstar` array with the solution values for the edges in the path.
	for (int i = 0; i < nnodes; i++) {
        int node1 = s->path[i];
        int node2 = s->path[i + 1];

		xstar[xpos(node1, node2, nnodes)] = 1.0;
	}
}