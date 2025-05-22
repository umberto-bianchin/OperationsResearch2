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

    extra_mileage_operation(inst, &s, nInserted, inserted);

    compute_solution_cost(inst, &s);

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

        double random = ((double)rand()/RAND_MAX);

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

/**
 * @brief Runs the genetic algorithm to solve the TSP instance within a time limit.
 *
 * Initializes a population of random tours, performs selection, crossover, repair,
 * and replacement until the time limit is reached. Updates the global best solution.
 *
 * @param inst Pointer to the TSP instance containing distances and parameters.
 * @param timelimit Maximum time (in seconds) to run the GA loop.
 */
void genetic_algorithm(instance *inst, double timelimit){
    int pop_size = inst->params[POPULATION_SIZE];
    int generation_size = inst->params[GENERATION_SIZE];

    solution *solutions = calloc(pop_size, sizeof(solution));
    solution *childs = calloc(generation_size, sizeof(solution));
    
    srand(inst->seed);

    // Generate random solutions
    generate_random_solutions(inst, solutions);

    // Apply genetic algorithm operations
    double elapsed_time = second() - inst->t_start;
    int iteration = 0;
    while(elapsed_time < timelimit){
        // Find the champion
        int best_idx = 0;          
        for(int i = 0; i < pop_size; i++){
            if(solutions[i].cost < solutions[best_idx].cost){
                best_idx = i;
            }
        }

        double best_cost = solutions[best_idx].cost;

        printf("Generation %d, best cost is %lf\n", iteration, solutions[best_idx].cost);
        printf("Generation %d, generating childs\n", iteration);

        // Create all the childs
        int generatedChilds = generate_childs(inst, solutions, childs, timelimit, best_cost);

        printf("Generation %d, generated %d childs\n", iteration, generatedChilds);

        //two_opt(inst, &inst->best_solution, inst->time_limit);

        printf("Generation %d, killing population members\n", iteration);

        // Calculate the total cost of the solutions
        double total_cost = 0;
        for(int i = 0; i < pop_size; i++){
            total_cost += solutions[i].cost;
        }

        // Kill the population (remove GENERATION_SIZE solutions), with a roulette-wheel elimination
        kill_members(inst, solutions, childs, total_cost, best_idx, generatedChilds);

        elapsed_time = second() - inst->t_start;
        iteration++;
    }

    int best_idx = 0;          
    for(int i = 0; i < pop_size; i++){
        if(solutions[i].cost < solutions[best_idx].cost){
            best_idx = i;
        }
    }

    printf("Final cost founded after %d generations is: %lf\n", iteration, solutions[best_idx].cost);
    update_best_solution(inst, &solutions[best_idx]);

    for(int i = 0; i < pop_size; i++){
        free_route(&solutions[i]);
    }
    for(int i = 0; i < generation_size; i++){
        free_route(&childs[i]);
    }

    free(childs);
    free(solutions);  
}

/**
 * @brief Inserts remaining nodes into a partial tour by minimal extra mileage.
 *
 * Starting from a tour prefix of length nInserted, repeatedly finds the best
 * insertion position for each unvisited node that minimizes the cost increase.
 *
 * @param inst Pointer to the TSP instance (contains cost matrix).
 * @param s Pointer to the solution being repaired (path[] must have prefix of length nInserted).
 * @param nInserted Number of nodes currently in the tour prefix.
 * @param inserted Boolean array of size nnodes marking which nodes are already in the prefix.
 */
void extra_mileage_operation(instance *inst, solution *s, int nInserted, int *inserted){
    int nodes = inst->nnodes;
    int a, b;
    int node_a = -1, node_b = -1, node_h = -1;
    double bestDelta, currentDelta;
    double distance, maxDist = -1;

    while(nInserted < nodes){ 
        bestDelta = INF_COST;
        node_a = -1; node_b = -1; node_h = -1;
        
        for(int pos = 0; pos < nInserted; pos++) {
            a = s->path[pos];
            b = s->path[pos + 1];
            
            for(int h = 0; h < nodes; h++){
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
            if(s->path[insertPos] == node_a && s->path[insertPos+1] == node_b)
                break;
        }
        
        for(int i = nInserted + 1; i > insertPos+1; i--){
            s->path[i] = s->path[i - 1];
        }
        
        s->path[insertPos+1] = node_h;
        inserted[node_h] = 1;
        nInserted++;
    }

    s->path[nodes] = s->path[0];
}

/**
 * @brief Removes duplicate nodes from a tour, preserving first occurrences.
 *
 * Takes the current path (possibly with repeats) and compacts it to a unique prefix.
 *
 * @param inst Pointer to the TSP instance (for nnodes).
 * @param s Pointer to the solution whose path will be compacted.
 * @return The number of unique nodes copied to the new path prefix.
 */
int short_cut(instance *inst, solution *s){
    // Arrays initialized at 0 through calloc
    int *new_path = calloc(inst->nnodes, sizeof(int));
    int *visited = calloc(inst->nnodes, sizeof(int));
    int nInserted = 0;
    for(int i = 0; i < inst->nnodes; i++){
        if(visited[s->path[i]] == 0){
            new_path[nInserted++] = s->path[i];
            visited[s->path[i]] = 1;
        }
    }

    memcpy(s->path, new_path, nInserted * sizeof(int));

    free(visited);
    free(new_path);
    return nInserted;
}

/**
 * @brief Generates an initial population of random tours for the TSP.
 *
 * Allocates and fills each solution.path with a random permutation of nodes
 * (plus return to the start), then evaluates its cost.
 *
 * @param inst Pointer to the TSP instance (contains nnodes and seed).
 * @param solutions Pre-allocated array of solutions of size pop_size.
 */
void generate_random_solutions(instance *inst, solution *solutions){
    int pop_size = inst->params[POPULATION_SIZE];

    // Generate random solutions
    srand(inst->seed);
    for(int i = 0; i < pop_size; i++){
        allocate_route(&solutions[i], inst->nnodes);
        int *choosen = calloc(inst->nnodes, sizeof(int));  // array initialized to 0 through calloc

        for(int j = 0; j < inst->nnodes; j++){
            int node;
            do {
                node = rand() % inst->nnodes;
            } while(choosen[node]);
            choosen[node] = 1;
            solutions[i].path[j] = node;
        }
        solutions[i].path[inst->nnodes] = solutions[i].path[0]; // close tour

        compute_solution_cost(inst, &solutions[i]);
        check_solution(inst, &solutions[i]);

        free(choosen);
    }
}

/**
 * @brief Generates child solutions via crossover, repair, and local improvement.
 *
 * Performs one-point crossover between pairs of parents, applies shortcut and
 * extra-mileage repair, then conditionally runs 2-opt on the best 5% of children.
 *
 * @param inst Pointer to the TSP instance.
 * @param solutions Parent population array.
 * @param childs Pre-allocated array to receive generated children.
 * @param timelimit Remaining time budget (to stop early if needed).
 * @param best_cost Cost of the current best solution for thresholding.
 * @return Number of children actually generated (<= generation_size).
 */
int generate_childs(instance *inst, solution *solutions, solution *childs, double timelimit, double best_cost){
    int generation_size = inst->params[GENERATION_SIZE];
    int pop_size = inst->params[POPULATION_SIZE];
    int generatedChilds = 0;

    for(int i = 0; i < generation_size; i++){
        allocate_route(&childs[i], inst->nnodes);

        // Select two random parents
        int parent1 = rand() % pop_size;
        int parent2 = rand() % pop_size;
        while(parent1 == parent2){
            parent2 = rand() % pop_size;
        }

        const solution a = solutions[parent1];
        const solution b = solutions[parent2];
        
        // Crossover operation
        int crossover_index = rand() % inst->nnodes;
        for(int j = 0; j < inst->nnodes; j++){
            if(j < crossover_index){
                childs[i].path[j] = a.path[j];
            } else {
                childs[i].path[j] = b.path[j];
            }
        }
        childs[i].path[inst->nnodes] = childs[i].path[0]; // Close the tour

        // Short-cut operation to remove duplicates
        int nInserted = short_cut(inst, &childs[i]);
        
        // Create the already visited nodes array
        int *inserted = calloc(inst->nnodes, sizeof(int));
        for(int k = 0; k <= nInserted; k++){
            inserted[childs[i].path[k]] = 1;
        }

        // Apply Extra mileage operation to get a valid solution
        extra_mileage_operation(inst, &childs[i], nInserted, inserted);

        compute_solution_cost(inst, &childs[i]);

        // Apply two opt to the childs that are at max 10% worse than the champion
        if (childs[i].cost < 1.05 * best_cost) {
            two_opt(inst, &childs[i], timelimit / 10);
        }

        check_solution(inst, &childs[i]);

        free(inserted);

        generatedChilds++;

        if(second() - inst->t_start > timelimit){
            break;
        }
    }

    return generatedChilds;
}

/**
 * @brief Eliminates and replaces population members via roulette-wheel.
 *
 * Removes `generatedChild` solutions from the population (except the champion)
 * based on cost-proportional probabilities, then copies in the corresponding
 * child solutions.
 *
 * @param inst Pointer to the TSP instance.
 * @param solutions Current population array.
 * @param childs Array of candidate child solutions.
 * @param total_cost Sum of all population solution costs.
 * @param best_idx Index of the champion (never eliminated).
 * @param generatedChild Number of valid children to insert.
 */
void kill_members(instance *inst, solution *solutions, solution *childs, double total_cost, int best_idx, int generatedChild){
    int generation_size = generatedChild;
    int pop_size = inst->params[POPULATION_SIZE];

    for(int i = 0; i < generation_size; i++){
        double r = ((double)rand() / RAND_MAX) * total_cost;
        int remove_idx = -1;

        // Walk through population costs subtracting until threshold reached
        for(int j = 0; j < pop_size; j++){
            if(j == best_idx) continue;
            r -= solutions[j].cost;
            if(r <= 0){
                remove_idx = j;
                break;
            }
        }

        // If no index found, remove the worst member
        if(remove_idx < 0){
            double worst = -1.0;
            for(int i = 0; i < pop_size; i++){
                if(i == best_idx) continue;
                if(solutions[i].cost > worst){
                    worst = solutions[i].cost;
                    remove_idx = i;
                }
            }
        }

        //free_route(&solutions[remove_idx]);
        memcpy(solutions[remove_idx].path, childs[i].path, (inst->nnodes + 1) * sizeof(int));
        solutions[remove_idx].cost = childs[i].cost;
    }
}