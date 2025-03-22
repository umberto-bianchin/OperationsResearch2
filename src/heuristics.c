#include <heuristics.h>
#include <utils.h>

/**
 * @brief
 * Compute the solution with the nearest neighbour heuristic algorithm given a starting node
 * @param inst the tsp instance
 * @param start_node the node used to start the algorithm
 */
void nearest_neighbour(instance *inst, int start_node){
    int nodes = inst->nnodes;
    int *visited = (int *) calloc(nodes, sizeof(int));

    inst->solution[0] = start_node;
    inst->solution[nodes] = start_node;

    visited[start_node] = 1;
    
    for(int i = 1; i < nodes; i++){
        int last_selected_node = inst->solution[i-1];
        int nearest_node = -1;
        double min_cost = INF_COST;

        for(int j = 0; j < nodes; j++){
            if(visited[j] == 0 && inst->costs[last_selected_node * nodes + j] < min_cost){
                min_cost = inst->costs[last_selected_node * nodes + j];
                nearest_node = j;
            }
        }

        if(nearest_node != -1){
            inst->solution[i] = nearest_node; 
            visited[nearest_node] = 1;
        }
        
    }
    
    free(visited);

    compute_solution_cost(inst);
    check_solution(inst);
}

/**
 * @brief
 * Compute the solution with the nearest neighbour heuristic algorithm analyzing all possible starting nodes
 * with respect to the time limit
 * @param inst the tsp instance
 */
void multi_start_nearest_neighbours(instance *inst){
    for(int i = 0; i < inst->nnodes; i++){
        nearest_neighbour(inst, i);
        two_opt(inst);
        check_solution(inst);

        double t2 = second();

        if(t2 - inst->t_start > inst->time_limit){
            if(VERBOSE>=ERROR){
                print_error("Exceded time limit while computing multi_start_nearest_neighbours, exiting the loop.\n", false);
                break;
            }
        }
    }
}

/**
 * @brief
 * Compute the solution with the variable veighbourhood algorithm, starting with a nearest neighbourhood solution
 * starting from a random node
 * @param inst the tsp instance
 */
void variable_neighbourhood(instance *inst){
    int iterations_without_improvement = 0;
    
    nearest_neighbour(inst, rand() % inst->nnodes);
    two_opt(inst);
    
    while (second() - inst->t_start < inst->time_limit &&
    iterations_without_improvement < MAX_NO_IMPROVEMENT) {
        
        double old_cost = inst->best_cost;
	    memcpy(inst->solution, inst->best_solution, (inst->nnodes + 1) * sizeof(int));

        for (int i = 0; i < KICK; i++) {
            three_opt(inst);
        }

        compute_solution_cost(inst);
        two_opt(inst);

        if(old_cost < inst->best_cost)
            iterations_without_improvement = 0;
        else
            iterations_without_improvement++;
    }
}

/**
 * @brief
 * Compute the solution with the extra mileage heuristic algorithm starting from the most distant points
 * @param inst the tsp instance
 */
void extra_mileage(instance *inst){
    double elapsed_time = second() - inst->t_start;

    int nodes = inst->nnodes;
    int i, j, a, b, h;
    int node_a = -1, node_b = -1, node_h = -1;
    int nInserted = 0; 
    double bestDelta, currentDelta;
    double distance, maxDist = -1;

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

    inst->solution[0] = node_a;
    inst->solution[1] = node_b;
    inserted[node_a] = 1;
    inserted[node_b] = 1;
    nInserted = 2;

    while(nInserted < nodes && elapsed_time < inst->time_limit){ 
        bestDelta = INF_COST;
        node_a = -1; node_b = -1; node_h = -1;
        
        for(int pos = 0; pos < nInserted; pos++) {
            a = inst->solution[pos];
            b = inst->solution[pos+1];
            
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
            if(inst->solution[insertPos] == node_a && inst->solution[insertPos+1] == node_b)
                break;
        }
        
        for(i = nInserted + 1; i > insertPos+1; i--){
            inst->solution[i] = inst->solution[i - 1];
        }
        
        inst->solution[insertPos+1] = node_h;
        inserted[node_h] = 1;
        nInserted++;

        elapsed_time = second() - inst->t_start;
    }

    inst->solution[inst->nnodes] = inst->solution[0];

    compute_solution_cost(inst);
    check_solution(inst);

    free(inserted);
}

/**
 * @brief
 * Greedy Randomized Adaptive Search Path algorithm
 * @param inst the tsp instance
 * @param start_node the starting node
 */
void grasp(instance *inst, int start_node) {
    int nodes = inst->nnodes;
    int *visited = (int *) calloc(nodes, sizeof(int));
    int *nearest_node = (int *) calloc(MIN_COSTS, sizeof(int));
    double *min_cost = (double *) calloc(MIN_COSTS, sizeof(double));

    inst->solution[0] = start_node;
    inst->solution[nodes] = start_node;

    visited[start_node] = 1;
    
    for(int i = 1; i < nodes; i++){
        int last_selected_node = inst->solution[i-1];

        for(int l = 0; l < MIN_COSTS; l++){
            min_cost[l] = INF_COST;
            nearest_node[l] = -1;
        }

        for(int j = 0; j < nodes; j++){
            if(visited[j] == 1)
                continue;
            
            double cost = inst->costs[last_selected_node * nodes + j];

            for(int k = 0; k < MIN_COSTS; k++){
                if(cost < min_cost[k]){
                    for(int l = MIN_COSTS - 1; l > k; l--){
                        min_cost[l] = min_cost[l-1];
                        nearest_node[l] = nearest_node[l-1];
                    }
                    
                    min_cost[k] = cost;
                    nearest_node[k] = j;
                    break;
                }
            }
        }

        double r = (double) rand() / (double) RAND_MAX;

        if(r <= ALPHA){
            int valid_count = 0;

            for(int k = 0; k < MIN_COSTS; k++){
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

                inst->solution[i] = best_node;
                visited[best_node] = 1;
            } else{
                int random_index = rand() % valid_count;
                inst->solution[i] = nearest_node[random_index];
                visited[nearest_node[random_index]] = 1;
            }
        } else {
            inst->solution[i] = nearest_node[0];
            visited[nearest_node[0]] = 1;
        }
    }
    
    free(visited);

    compute_solution_cost(inst);
    check_solution(inst);
}

/**
 * @brief
 * Compute the solution with the Greedy Randomized Adaptive Search Path + Two Opt local search
 * algorithm analyzing all possible starting nodes with respect to the time limit
 * @param inst the tsp instance
 */
void multi_start_grasp(instance *inst) {    
    double t2;

    for(int i = 0; i < inst->nnodes; i++){
        grasp(inst, i);
        two_opt(inst);
        check_solution(inst);

        t2 = second();
        
        if(t2 - inst->t_start > inst->time_limit){
            if(VERBOSE>=ERROR){
                print_error("Exceded time limit while computing multi_start_grasp, exiting the loop.\n", false);
                break;
            }
        }
    }
}

/**
 * @brief
 * Compute the solution with the TABU search algorithm
 * @param inst the tsp instance
 */
void tabu(instance *inst){
    int nodes = inst->nnodes;

    // Initialize all nodes as non-tabu
    int **tabuList = malloc(nodes * sizeof(int *));
    for (int i = 0; i < nodes; i++) {
        tabuList[i] = calloc(nodes, sizeof(int));
    }

    int currentTenure = MIN_TENURE;

    nearest_neighbour(inst, rand() % inst->nnodes);
    int iter = 0;

	while (((second() - inst->t_start) < inst->time_limit)){	
		double min_delta = INF_COST;
		int swap_i = -1, swap_j = -1;
		
		for (int i = 1; i < nodes; i++) {
			for (int j = i + 2; j < nodes; j++) {

                if (iter < tabuList[inst->solution[i]][inst->solution[j]])
                    continue;

				double current_delta = calculate_delta(i, j, inst);

                if(current_delta < min_delta){
					min_delta = current_delta;
					swap_i = i;
					swap_j = j;
				}
			}
		}

        if(swap_i == -1 || swap_j == -1){
            if(VERBOSE >= ERROR)
                print_error("No possible swap found, exiting the loop\n", false);

            break;
        }
        
        tabuList[inst->solution[swap_i]][inst->solution[swap_j]] = iter + currentTenure;

        if(VERBOSE >= DEBUG)
            printf("Swapping node %d with node %d\n", swap_i, swap_j);

        reverse_segment(swap_i + 1, swap_j, inst);
        compute_solution_cost(inst);
        check_solution(inst);

        double elapsed_time = second() - inst->t_start;

        if(elapsed_time > inst->time_limit){
            if(VERBOSE>=ERROR)
                print_error("Exceded time limit while computing tabu search, exiting the loop\n", false);

            break;
        }

        iter++;
        currentTenure += TENURE_STEP;
        if (currentTenure > MAX_TENURE)
            currentTenure = MIN_TENURE;
    }

    for (int i = 0; i < nodes; i++) {
        free(tabuList[i]);
    }
    free(tabuList);
}