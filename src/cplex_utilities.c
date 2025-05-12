#define _CRT_SECURE_NO_WARNINGS
#include <cplex_utilities.h>

/**
 * @brief Compute index of binary variable x(i,j) in lexicographic ordering for i<j.
 *
 * Ordering is based on flattening upper triangle of adjacency matrix.
 * @param i       First node index (0-based).
 * @param j       Second node index (0-based).
 * @param nnodes  Number of nodes in graph.
 * @return        Position of x(i,j) in CPLEX variable array.
 */
int xpos(int i, int j, int nnodes){
	if (i == j) print_error(" i == j in xpos" );
	if (i > j) return xpos(j,i,nnodes);
	int pos = i * nnodes + j - ((i + 1) * (i + 2)) / 2;

	return pos;
}

/**
 * @brief Build CPLEX model: add variables and degree constraints for TSP.
 *
 * @param inst  Pointer to TSP instance with coords and nnodes.
 * @param env   CPLEX environment pointer.
 * @param lp    CPLEX problem pointer.
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
 * @brief Solve TSP via CPLEX Benders or Branch-and-Cut or Concorde callback.
 *
 * Builds model, sets parameters, optionally warms up, and iterates
 * adding SECs until one component. Uses patching heuristic if time.
 * @param inst  TSP instance with time_limit, algorithm flag.
 * @return      0 on success (exits on error internally).
 */
int TSPopt(instance *inst){  
	// Open CPLEX model
	int error = 0;
	CPXENVptr env = CPXopenCPLEX(&error);
	if (error) print_error("CPXopenCPLEX() error");
	CPXLPptr lp = CPXcreateprob(env, &error, "TSP model version 1"); 
	if (error) print_error("CPXcreateprob() error");

	build_model(inst ,env, lp);
	
	// Cplex's parameter setting
	CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_OFF);
	CPXsetdblparam(env, CPX_PARAM_EPGAP, 1e-9);
	if(VERBOSE >= DEBUG) CPXsetintparam(env, CPX_PARAM_SCRIND, CPX_ON); // Cplex output on screen
	
	int *comp = (int *) calloc(inst->nnodes, sizeof(int));
	int *succ = (int *) calloc(inst->nnodes, sizeof(int));
	inst->ncols = CPXgetnumcols(env, lp);
	double *xstar = (double *) calloc(inst->ncols, sizeof(double));
	int ncomp = 9999;
	double objval = CPX_INFBOUND;
	int iter = 0;
	
	if(inst->params[WARMUP]){
		warmup_CPX_solution(inst, env, lp);
	}

	// register callback for branch-and-cut
	if(inst->algorithm == 'C'){
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

	// main loop: solve MIP, extract x*, separate subtours
	while(ncomp >= 2){
		double elapsed_time = second() - inst->t_start;
		double residual_time = inst->time_limit - elapsed_time;
		CPXsetdblparam(env, CPX_PARAM_TILIM, residual_time);
		
		error = CPXmipopt(env,lp);

		if (error){
			printf("CPX error code %d\n", error);
			print_error("CPXmipopt() error"); 
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

		elapsed_time = second() - inst->t_start;
		residual_time = inst->time_limit - elapsed_time;

		if(residual_time <= 0){
			printf("Exceded time limit while computing TSPopt(), exiting the loop\n");
			break;
		}
	}
	
	// Write the model in an appropriate file 
	if (VERBOSE >= DEBUG)
		CPXwriteprob(env, lp, "history/model.lp", NULL);

	// optional patching heuristic if subtours remain
	if(ncomp >= 2){
		if(VERBOSE >= INFO){
			printf("Entering Patching Heuristic method with ncomp %4d, time %5.2lf\n", ncomp, second() - inst->t_start);
			fflush(NULL);
		}

		patching_heuristic(inst, succ, comp, ncomp);

		if(VERBOSE >= INFO){
			printf("Exiting Patching Heuristic method at time %5.2lf\n", second() - inst->t_start);
			fflush(NULL);
		}
	}

	// recover best subtour-free solution
	solution s;
	allocate_route(&s, inst->nnodes);
	solution_from_CPX(inst, &s, succ);
	update_best_solution(inst, &s);

	// Free and close cplex model
	free_route(&s);
	free(xstar);
	free(succ);
	free(comp);
	CPXfreeprob(env, &lp);
	CPXcloseCPLEX(&env); 

	return 0;
}

/**
 * @brief Build succ[] and comp[] arrays from CPLEX solution xstar[].
 *
 * Identifies connected components (subtours) and assigns
 * succ[i] = next node in tour, comp[i] = component ID.
 */
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
 * Dispatch CPLEX callback based on context: candidate cuts or relaxation cuts.
 *
 * Determines whether to call sec_callback (lazy constraints) or
 * concorde_callback (user cuts via Concorde routines).
 *
 * @param context    CPLEX callback context pointer.
 * @param contextid  Identifier of callback context (CANDIDATE or RELAXATION).
 * @param userhandle User data (cast to instance*).
 * @return 0 on success; non-zero aborts solve.
 */
int CPXPUBLIC cpx_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle){ 
	instance* inst = (instance*) userhandle;
	if(contextid == CPX_CALLBACKCONTEXT_CANDIDATE) return sec_callback(context, contextid, inst); 
	if(contextid == CPX_CALLBACKCONTEXT_RELAXATION) return concorde_callback(context, contextid, inst);
	print_error("Contextid unknownn in cpx_callback()"); 
	return 1; 	
}

/**
 * @brief
 * Candidate callback for branch-and-cut: checks integer candidate solution,
 * rejects if subtours exist by adding lazy SECs, and optionally posts heuristic.
 *
 * @param context    CPLEX callback context pointer.
 * @param contextid  Must be CPX_CALLBACKCONTEXT_CANDIDATE.
 * @param userhandle User data (cast to instance*).
 * @return 0 on success.
 */
int CPXPUBLIC sec_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle){ 
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

	// If more than one component, add lazy SECs
	if(ncomp >= 2){
		add_sec(inst, NULL, NULL, context, contextid, comp, ncomp, true);
	}
	
	double elapsed_time = second() - inst->t_start;
	add_solution(&inst->history_best_costs, objval, elapsed_time);

	// Optionally post a heuristic solution at certain nodes
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
			post_CPX_solution(inst, context, succ, comp, ncomp);
		}
	}
	
	free(succ);
	free(comp);
	free(xstar);
	return 0; 
}

/**
 * @brief
 * Relaxation callback using Concorde's mincut routines to generate user cuts (SECs).
 * Only invoked at the root node to limit overhead.
 *
 * @param context    CPLEX callback context pointer.
 * @param contextid  Must be CPX_CALLBACKCONTEXT_RELAXATION.
 * @param userhandle User data (cast to instance*).
 * @return 0 on success.
 */
int CPXPUBLIC concorde_callback(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *userhandle){
	instance* inst = (instance*) userhandle;  
	int current_node = -1;

	if(CPXcallbackgetinfoint(context, CPXCALLBACKINFO_NODECOUNT, &current_node)){
		print_error("CPXcallbackgetinfopoint error");
	}

	// Only at first relaxation (root) to save time
	if (current_node % inst->nnodes != 0){
		return 0;
	}

	int* elist =(int *) calloc(2 * inst->ncols, sizeof(int)); // elist contains each pair of vertex such as (1,2), (1,3), (1,4), (2, 3), (2,4), (3,4) so in list becomes: 1,2,1,3,1,4,2,3,2,4,3,4
	int *comps = NULL;
    int *compscount = NULL;
    int num_edges = 0, k = 0, ncomp = 0;
	double* xstar = (double*) malloc(inst->ncols * sizeof(double));  
	double objval = CPX_INFBOUND;

	for (int i = 0; i < inst->nnodes; i++){
		for (int j = i + 1; j < inst->nnodes; j++){
			elist[k++] = i;
			elist[k++] = j;
			num_edges++;
		}
	}

	if (CPXcallbackgetrelaxationpoint(context, xstar, 0, inst->ncols-1, &objval)){
		print_error("CPXcallbackgetrelaxationpoint error");
	}

	// Checking whether or not the graph is connected with the fractional solution.
	if (CCcut_connect_components(inst->nnodes, num_edges, elist, xstar, &ncomp, &compscount, &comps)){
		print_error("CCcut_connect_components() error ");
	}

	if (ncomp == 1){
		// Graph connected: generate violated cuts via Concorde
		double cutOff = 1.9; //required a violation of at least 0.1 in the sec in the cut form
		relaxation_callback_params params = { .context = context, .inst = inst, .xstar = xstar };

		if (CCcut_violated_cuts(inst->nnodes, num_edges, elist, xstar, cutOff, violated_cuts_callback, &params)){
			print_error("%d, CCcut_violated_cuts() error ");
		}
	} else {
		// Multiple subtours: add SECs via lazy mechanism
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
 * Concorde user-cut callback: transforms detected violated SECs into CPLEX cuts.
 *
 * @param cutval      Minimum required violation threshold.
 * @param n           Number of nodes in the cut.
 * @param members     List of node indices in the cut.
 * @param userhandle  relaxation_callback_params pointer.
 * @return            CPXcallbackaddusercuts() status.
 */
int violated_cuts_callback(double cutval, int n, int *members, void *userhandle){
    relaxation_callback_params *params = (relaxation_callback_params *)userhandle;
    instance *inst = params->inst;
    CPXCALLBACKCONTEXTptr context = params->context;

    double rhs = n - 1;
    char sense = 'L';
    int matbeg = 0;
	int purgeable = CPX_USECUT_FILTER;
    int local = 0;

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
    
    int status = CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &matbeg, index, value, &purgeable, &local);

	/*double violation = cut_violation(nnz, rhs, sense, matbeg, index, value, params->xstar);
	double expected_violation = (2.0 - cutval) / 2.0;
	
	if(VERBOSE >= DEBUG){
		printf("Generated a violated cut: violation %lf, expected violation %lf\n", violation, expected_violation);
		fflush(NULL);
	}

	if(fabs(violation - expected_violation) > 0.0001){
		print_error("Cut not violated\n");
	}*/

    free(index);
    free(value);
    return status;
}

/**
 * @brief
 * Add or reject Subtour Elimination Constraints (SECs) based on subtour components.
 *
 * When in callback context, rejects candidate or posts user cuts;
 * otherwise adds rows directly to the model.
 *
 * @param inst       TSP instance containing problem data.
 * @param env        CPLEX environment (NULL in callback).
 * @param lp         CPLEX problem pointer (NULL in callback).
 * @param context    Callback context pointer (NULL outside callback).
 * @param contextid  Callback context ID (RELAXATION or CANDIDATE).
 * @param comp       Array of component IDs for each node (length nnodes).
 * @param ncomp      Number of connected components (>1 indicates subtours).
 * @param callback   True if inside callback, false if modeling phase.
 */
void add_sec(instance *inst, CPXENVptr env, CPXLPptr lp, CPXCALLBACKCONTEXTptr context, CPXLONG contextid, int *comp, int ncomp, bool callback){
	int izero = 0; 
	char **cname = (char **) calloc(1, sizeof(char *));
	int max_edges = inst->nnodes * (inst->nnodes - 1) / 2;

	int *index = (int *) calloc(max_edges, sizeof(int));
	double *value = (double *) calloc(max_edges, sizeof(double));
	cname[0] = (char *) calloc(100, sizeof(char));

	// Iterate each component to build its SEC
	for(int k = 1; k <= ncomp; k++){
		double rhs = -1;
		char sense = 'L';
		int nnz = 0;

		// Count members and edges within component k
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

		// Add or reject based on context
		if(callback){
			if(contextid == CPX_CALLBACKCONTEXT_CANDIDATE){
				// Lazy cut rejection
				if(CPXcallbackrejectcandidate(context, 1, nnz, &rhs, &sense, &izero, index, value)){
					print_error("CPXcallbackrejectcandidate() error");
				}
			} else if(contextid == CPX_CALLBACKCONTEXT_RELAXATION){
				// User cut on relaxation
				int purgeable = CPX_USECUT_FILTER;
				int local = 0;

				if(CPXcallbackaddusercuts(context, 1, nnz, &rhs, &sense, &izero, index, value, &purgeable, &local)){
					print_error("CPXcallbackaddusercuts() error");
				}
			}
		} else{
			// Modeling phase: add SEC row directly
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
 * Convert CPLEX successor array into tour path and compute its cost.
 *
 * @param inst  TSP instance (contains cost matrix, nnodes).
 * @param s     Solution struct with path[] allocated length nnodes+1.
 * @param succ  Successor array: succ[i] = next node in tour.
 */
void solution_from_CPX(instance *inst, solution *s, int *succ) {
    int current = 0; 
    for (int i = 0; i <= inst->nnodes; i++) {
		s->path[i] = current;
        current = succ[current];
    }

	double total_cost = 0.0;
	for (int i = 0; i < inst->nnodes; i++)
		total_cost += inst->costs[s->path[i] * inst->nnodes + s->path[i +1]];
	

	s->cost = total_cost;
}

/**
 * @brief
 * Heuristic patching to merge subtours into single tour by minimal reconnection.
 *
 * Edges between components are swapped based on minimal delta cost.
 *
 * @param inst   TSP instance with cost matrix.
 * @param succ   Successor array modified in-place.
 * @param comp   Component ID array for each node (updated in-place).
 * @param ncomp  Number of components before merging.
 */
void patching_heuristic(instance *inst, int *succ, int *comp, int ncomp){
    if (ncomp < 2){
        return;
    }

	// For each pair of distinct components
	for (int c1 = 1; c1 <= ncomp; c1++){
		for (int c2 = c1 + 1; c2 <= ncomp; c2++){

			// skip if second comp is empty
			bool valid_c2 = false;
            for (int k = 0; k < inst->nnodes; k++) {
                if (comp[k] == c2) {
                    valid_c2 = true;
                    break;
                }
            }
            if (!valid_c2) continue;

			// Collect nodes in each comp
			int *nodes_c1 = (int *) calloc(inst->nnodes, sizeof(int));
			int *nodes_c2 = (int *) calloc(inst->nnodes, sizeof(int));
			int size_c1 = 0, size_c2 = 0;

			for (int i = 0; i < inst->nnodes; i++){
				if (comp[i] == c1){
					nodes_c1[size_c1++] = i;
				}
				else if(comp[i] == c2) {
					nodes_c2[size_c2++] = i;
				}
			}

			// Find best reconnection with minimal delta
			int best_i = -1, succ_i = -1, best_j = -1, succ_j = -1;
			double min_delta = INF_COST;
			bool orientation;
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
				reverse_cycle(best_j, succ);

				succ[best_i] = best_j;
				succ[succ_j] = succ_i;
			}
			else{
				succ[best_i] = succ_j;
				succ[best_j] = succ_i;
			}

			// Merge labels
			for(int j = 0; j < size_c2; j++)
				comp[nodes_c2[j]] = c1;

			
			free(nodes_c2);
			free(nodes_c1);
		}
	}
}

/**
 * @brief
 * Compute cost delta for merging edges between two components.
 *
 * @param inst    TSP instance with cost matrix.
 * @param i1,j1   Edge in comp1 and comp2.
 * @param i2,j2   Successor edges in comps.
 * @param option  true=swap cross-edges (i1->j1 & i2->j2),
 *                false=swap alt-edges (i1->j2 & i2->j1).
 * @return        Difference in cost (new - old).
 */
double delta_cost(instance *inst, int i1, int j1, int i2, int j2, bool option){
	if(option)
		// reconnect as (i1->j1) + (i2->j2)
    	return (inst->costs[i1 * inst->nnodes + j1] + inst->costs[i2 * inst->nnodes + j2]) - (inst->costs[i1 * inst->nnodes + i2] + inst->costs[j1 * inst->nnodes + j2]);
	else
		// reconnect as (i1->j2) + (i2->j1)
		return (inst->costs[i1 * inst->nnodes + j2] + inst->costs[i2 * inst->nnodes + j1]) - (inst->costs[i1 * inst->nnodes + i2] + inst->costs[j1 * inst->nnodes + j2]);	
}


/**
 * @brief
 * Reverse a cycle in the successor representation, starting at `start`.
 *
 * @param start  Node to start cycle reversal.
 * @param succ   Successor array modified in-place.
 */
void reverse_cycle(int start, int *succ){
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
 * Generate initial heuristic solution and inject it as CPLEX MIP start.
 *
 * Uses nearest neighbour + 2-opt for a fraction of time_limit.
 *
 * @param inst  TSP instance with time_limit, best_solution set.
 * @param env   CPLEX environment pointer.
 * @param lp    CPLEX problem pointer.
 */
void warmup_CPX_solution(instance *inst, CPXENVptr env, CPXLPptr lp) {
	printf("Initializing solution with heuristics\n");

	double initialization_timelimit = inst->time_limit/10.0;
	nearest_neighbour(inst, rand() % inst->nnodes);     //need to initialize and optimize the first solution
	solution s;
	copy_solution(&s, &inst->best_solution, inst->nnodes);
	two_opt(inst, &s, initialization_timelimit);


    // Check if we have a valid solution to use as warm start
    if (s.path == NULL) {
        print_error("No solution available in s.path\n");
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
    printf("Initial solution found with cost %lf after %3.2lf seconds\n", s.cost, (second() - inst->t_start));
    free(value);
    free(index);
	free_route(&s);
}

/**
 * @brief
 * Post a heuristic solution from within a callback if it's better than incumbent.
 *
 * Extracts current solution, applies patching & local search, posts via callback API.
 *
 * @param inst     TSP instance.
 * @param context  Callback context pointer.
 * @param succ     Successor array from CPLEX solution.
 * @param comp     Component ID array.
 * @param ncomp    Number of components.
 */
void post_CPX_solution(instance *inst, CPXCALLBACKCONTEXTptr context, int *succ, int *comp, int ncomp){
	double elapsed_time = second() - inst->t_start;
	double residual_time = (inst->time_limit) - elapsed_time;

	patching_heuristic(inst, succ, comp, ncomp);

	solution s;
	allocate_route(&s, inst->nnodes);
	solution_from_CPX(inst, &s, succ);

	two_opt(inst, &s, residual_time);

	double incumbent = CPX_INFBOUND;
	double x_dummy;
	int hasIncumbent = 0;

	int error = CPXcallbackgetinfoint(context, CPXCALLBACKINFO_FEASIBLE, &hasIncumbent); 

	if(error){
		print_error("CPXcallbackgetinfoint() error");
	}
	if (!hasIncumbent) {
		return;
	}

	error = CPXcallbackgetincumbent(context, &x_dummy, 0, 0, &incumbent);
	if(error){
		char errmsg[CPXMESSAGEBUFSIZE];
		CPXgeterrorstring(NULL, error, errmsg);
		fprintf(stderr, "CPXcallbackgetincumbent Error %d: %s\n", error, errmsg);
		print_error("CPXcallbackgetincumbent() error");
	}

	int max_edges = inst->ncols;
	int *index = (int *) calloc(max_edges, sizeof(int));
	double *xstar = (double *) calloc(max_edges, sizeof(double));

	solution_to_CPX(&s, inst->nnodes, index, xstar);

	error = 0;
	if(s.cost < incumbent - EPS_COST){
		error = CPXcallbackpostheursoln(context, max_edges, index, xstar, s.cost, CPXCALLBACKSOLUTION_NOCHECK);
		if(VERBOSE >= INFO && !error){
			printf("Posted solution of cost %lf (my incubment was %lf)\n", s.cost, incumbent);
			fflush(NULL);
		} 
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
 * This function translates the solution stored in `s.path` into the format
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

/**
 * @brief
 * Compute violation amount of a candidate cut row.
 *
 * @param nnz    Number of nonzeros in cut.
 * @param rhs    Right-hand-side of the cut.
 * @param sense  'L','G','E' indicating <=,>=,=.
 * @param matbeg Always 0 for single row.
 * @param index  Array of variable positions.
 * @param value  Array of coefficients.
 * @param xstar  Current fractional solution values.
 * @return       Positive violation magnitude; zero if satisfied.
 */
double cut_violation(int nnz, double rhs, char sense, int matbeg, int *index, double *value, double *xstar){
	double cut_val = 0.0;
	for(int i = 0; i < nnz; i++){
		cut_val += value[i] * xstar[index[i]];
	}

	if(sense == 'L'){
		return fmax(cut_val - rhs, 0.0); //if the cut is not violated, this value will be negative
	} else if(sense == 'G'){
		return fmax(rhs - cut_val, 0.0);
	} else {
		return fabs(rhs - cut_val);
	}

	return 0;
}