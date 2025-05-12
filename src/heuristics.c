#include <heuristics.h>
#include <utils.h>

/**
 * @brief Build a tour using the nearest neighbour heuristic from a given start.
 *
 * @param inst        Pointer to TSP instance with cost matrix and nnodes set.
 * @param start_node  Index of node to start the greedy construction.
 */
void nearest_neighbour(instance *inst, int start_node){
    int nodes = inst->nnodes;
    // visited flags initialised to 0 via calloc
    int *visited = (int *) calloc(nodes, sizeof(int));

    solution s;
    // works on a local solution for thread safety
    copy_solution(&s, &inst->best_solution, inst->nnodes);

    s.path[0] = start_node;
    s.path[nodes] = start_node;

    visited[start_node] = 1;
    
    for(int i = 1; i < nodes; i++){
        int last_selected_node = s.path[i-1];
        int nearest_node = -1;
        double min_cost = INF_COST;

        for(int j = 0; j < nodes; j++){
            if(visited[j] == 0 && inst->costs[last_selected_node * nodes + j] < min_cost){
                min_cost = inst->costs[last_selected_node * nodes + j];
                nearest_node = j;
            }
        }

        if(nearest_node != -1){
            s.path[i] = nearest_node; 
            visited[nearest_node] = 1;
        }
    }
    
    
    compute_solution_cost(inst, &s);
    check_solution(inst, &s);
    update_best_solution(inst, &s);
    
    free(visited);
    free_route(&s);
}

/**
 * @brief Repeated nearest neighbour runs from each start within time budget.
 *
 * @param inst       Pointer to TSP instance with nnodes and t_start set.
 * @param timelimit  Maximum allowed CPU time for all restarts.
 */
void multi_start_nearest_neighbours(instance *inst, double timelimit){
    for(int i = 0; i < inst->nnodes; i++){
        nearest_neighbour(inst, i);

        // local copy to apply 2-opt without clobbering best
        solution s;
        allocate_route(&s, inst->nnodes);
        copy_solution(&s, &inst->best_solution, inst->nnodes);

        two_opt(inst, &s, timelimit);
        check_solution(inst, &(inst->best_solution));
        update_best_solution(inst, &s);

        if(second() - inst->t_start > timelimit){
            if(VERBOSE>=ERROR)
                printf("Exceded time limit while computing multi_start_nearest_neighbours, exiting the loop.\n");
            
            break;
        }

        free_route(&s);
    }
}

/**
 * @brief Variable Neighbourhood Search: apply perturbation+k-opt until no improvement.
 *
 * @param inst       Pointer to TSP instance with parameters KICK, K_OPT.
 * @param timelimit  CPU time limit for the search.
 */
void variable_neighbourhood(instance *inst, double timelimit){
    int iterations_without_improvement = 0;
    solution s;
    allocate_route(&s, inst->nnodes);
    copy_solution(&s, &inst->best_solution, inst->nnodes);
    
    while (second() - inst->t_start < timelimit &&
    iterations_without_improvement < MAX_NO_IMPROVEMENT) {

        double old_cost = inst->best_solution.cost;
         // reset to best before kicks
        memcpy(s.path, inst->best_solution.path, (inst->nnodes + 1) * sizeof(int));
        
        // apply KICK random k-opt moves
        for (int i = 0; i < inst->params[KICK]; i++) {
            if(inst->params[K_OPT] == 3){
                three_opt(inst, &s);
            }else if(inst->params[K_OPT] == 5){
                five_opt(inst, &s);
            }else {
                random_k_opt(inst, &s, inst->params[K_OPT]);
            }
        }

        compute_solution_cost(inst, &s);
        two_opt(inst, &s, timelimit);
        update_best_solution(inst, &s);

        // reset no-improv counter if improved
        if(old_cost < inst->best_solution.cost)
            iterations_without_improvement = 0;
        else
            iterations_without_improvement++;
    }

    free_route(&s);
}

/**
 * @brief Construct tour by repeatedly inserting farthest-increase node (extra mileage).
 *
 * @param inst Pointer to TSP instance; nnodes and costs must be set.
 */
void extra_mileage(instance *inst){
    double elapsed_time = second() - inst->t_start;
    int nodes = inst->nnodes;
    int i, j, a, b, h;
    int node_a = -1, node_b = -1, node_h = -1;
    int nInserted = 0; 
    double bestDelta, currentDelta;
    double distance, maxDist = -1;
    
    solution s;
    copy_solution(&s, &inst->best_solution, inst->nnodes);
    
    // track inserted nodes, zero-init
    int *inserted = (int *) calloc(nodes, sizeof(int));

    for(i = 0; i < nodes; i++){
        for(j = i + 1; j < nodes; j++){
            distance = inst->costs[i * nodes + j];
            if(distance > maxDist){
                maxDist = distance;
                node_a = i;
                node_b = j;
            }
        }
    }

    s.path[0] = node_a;
    s.path[1] = node_b;
    inserted[node_a] = 1;
    inserted[node_b] = 1;
    nInserted = 2;

    while(nInserted < nodes){ 
        bestDelta = INF_COST;
        node_a = -1; node_b = -1; node_h = -1;
        
        for(int pos = 0; pos < nInserted; pos++) {
            a = s.path[pos];
            b = s.path[pos + 1];
            
            for(h = 0; h < nodes; h++){
                if(inserted[h] == 0) {
                    currentDelta = inst->costs[a * nodes + h] + inst->costs[h * nodes + b] - inst->costs[a * nodes + b];
                    if(currentDelta < bestDelta){
                        bestDelta = currentDelta;
                        node_a = a;
                        node_b = b;
                        node_h = h;
                    }
                }
            }
        }
    
        if(node_h == -1){
            break;
        }

        int insertPos = 0;
        for(insertPos = 0; insertPos < nInserted; insertPos++){
            if(s.path[insertPos] == node_a && s.path[insertPos+1] == node_b)
                break;
        }
        
        for(i = nInserted + 1; i > insertPos+1; i--){
            s.path[i] = s.path[i - 1];
        }
        
        s.path[insertPos+1] = node_h;
        inserted[node_h] = 1;
        nInserted++;
    }

    s.path[inst->nnodes] = s.path[0];

    compute_solution_cost(inst, &s);
    update_best_solution(inst, &s);

    free_route(&s);
    free(inserted);
}

/**
 * @brief Greedy Randomized Adaptive Search Procedure (GRASP) constructive phase.
 *
 * @param inst        TSP instance with ALPHA and MIN_COSTS parameters.
 * @param start_node  Starting node index.
 */
void grasp(instance *inst, int start_node) {
    int nodes = inst->nnodes;
    int *visited = (int *) calloc(nodes, sizeof(int));
    int *nearest_node = (int *) calloc(inst->params[MIN_COSTS], sizeof(int));
    double *min_cost = (double *) calloc(inst->params[MIN_COSTS], sizeof(double));

    solution s;
    copy_solution(&s, &(inst->best_solution), inst->nnodes);

    s.path[0] = start_node;
    s.path[nodes] = start_node;

    visited[start_node] = 1;
    
    for(int i = 1; i < nodes; i++){
        int last_selected_node = s.path[i-1];

        for(int l = 0; l < inst->params[MIN_COSTS]; l++){
            min_cost[l] = INF_COST;
            nearest_node[l] = -1;
        }

        for(int j = 0; j < nodes; j++){
            if(visited[j] == 1)
                continue;
            
            double cost = inst->costs[last_selected_node * nodes + j];

            for(int k = 0; k < inst->params[MIN_COSTS]; k++){
                if(cost < min_cost[k]){
                    for(int l = inst->params[MIN_COSTS] - 1; l > k; l--){
                        min_cost[l] = min_cost[l-1];
                        nearest_node[l] = nearest_node[l-1];
                    }
                    
                    min_cost[k] = cost;
                    nearest_node[k] = j;
                    break;
                }
            }
        }

        double random = (double) rand() / (double) RAND_MAX;

        // alpha param must be divided by 100.0 since it is stored as integer 
        if(random <= inst->params[ALPHA] / 100.0){
            int valid_count = 0;

            for(int k = 0; k < inst->params[MIN_COSTS]; k++){
                if(nearest_node[k] != -1) {
                    valid_count++;
                } else {
                    break;
                }
            }

            if(valid_count == 0){
                double best_cost = INF_COST;
                int best_node = -1;

                for(int j = 0; j < nodes; j++){
                    if(!visited[j]) {
                        double c = inst->costs[last_selected_node*nodes + j];
                        if(c < best_cost) {
                            best_cost = c;
                            best_node = j;
                        }
                    }
                }

                s.path[i] = best_node;
                visited[best_node] = 1;
            } else{
                int random_index = rand() % valid_count;
                s.path[i] = nearest_node[random_index];
                visited[nearest_node[random_index]] = 1;
            }
        } else {
            s.path[i] = nearest_node[0];
            visited[nearest_node[0]] = 1;
        }
    }
    
    compute_solution_cost(inst, &s);
    check_solution(inst, &s);
    update_best_solution(inst, &s);

    free(visited);
    free_route(&s);
}

/**
 * @brief Multi-start GRASP: repeat GRASP + 2-opt across all starts within time.
 *
 * @param inst       TSP instance with time limit set.
 * @param timelimit  CPU time budget.
 */
void multi_start_grasp(instance *inst, double timelimit) {    
    for(int i = 0; i < inst->nnodes; i++){
        grasp(inst, i);
        
        solution s;
        allocate_route(&s, inst->nnodes);
        copy_solution(&s, &inst->best_solution, inst->nnodes);
        two_opt(inst, &s, timelimit);
        check_solution(inst, &inst->best_solution);
        update_best_solution(inst, &s);

        if(second() - inst->t_start > timelimit){
            if(VERBOSE>=ERROR){
                printf("Exceded time limit while computing multi_start_grasp, exiting the loop.\n");
                break;
            }
        }

        free_route(&s);
    }
}

/**
 * @brief Tabu Search: iteratively swap non-tabu edges with dynamic tenure.
 *
 * @param inst       TSP instance with TENURE and MAX_TENURE params.
 * @param timelimit  CPU time limit.
 */
void tabu(instance *inst, double timelimit){
    int nodes = inst->nnodes;
    solution s;
    allocate_route(&s, nodes);
    copy_solution(&s, &(inst->best_solution), nodes);
    
    // Tabu tenure matrix, zero-init via calloc
    int **tabuList = malloc(nodes * sizeof(int *));
    for (int i = 0; i < nodes; i++) {
        tabuList[i] = calloc(nodes, sizeof(int));
    }

    int currentTenure = inst->params[MIN_TENURE];
    int iter = 0;

	while (((second() - inst->t_start) < timelimit)){	
		double min_delta = INF_COST;
		int swap_i = -1, swap_j = -1;
		
		for (int i = 0; i < nodes; i++) {
			for (int j = i + 2; j < nodes; j++) {

                if (iter < tabuList[s.path[i]][s.path[j]])
                    continue;

				double current_delta = calculate_delta(i, j, inst, &s);

                if(current_delta < min_delta){
					min_delta = current_delta;
					swap_i = i;
					swap_j = j;
				}
			}
		}

        if(swap_i == -1 || swap_j == -1)
            print_error("No possible swap found, exiting the loop\n");
        
        // mark tabu for reversed edge
        tabuList[s.path[swap_i]][s.path[swap_j]] = iter + currentTenure;
        tabuList[s.path[swap_j]][s.path[swap_i]] = iter + currentTenure;

        if(VERBOSE >= DEBUG)
            printf("Swapping node %d with node %d\n", swap_i, swap_j);

        reverse_segment(swap_i, swap_j, &s);
        compute_solution_cost(inst, &s);
        check_solution(inst, &s);
        update_best_solution(inst, &s);

        iter++;
        currentTenure += inst->params[TENURE_STEP];
        // bound tenure
        if (currentTenure > inst->params[MAX_TENURE])
            currentTenure = inst->params[MIN_TENURE];
    }

    free_route(&s);
    for (int i = 0; i < nodes; i++) {
        free(tabuList[i]);
    }
    free(tabuList);
}