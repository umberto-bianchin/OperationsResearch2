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

    inst->solution_cost = compute_solution_cost(inst, inst->solution);
    check_solution(inst, false);
    two_opt(inst);
    check_solution(inst, false);
}

/**
 * @brief
 * Compute the solution with the nearest neighbour heuristic algorithm analyzing all possible starting nodes
 * with respect to the time limit
 * @param inst the tsp instance
 */
void multi_start_nearest_neighbours(instance *inst){
    double t1 = second();
    
    for(int i = 0; i < inst->nnodes; i++){
        nearest_neighbour(inst, i);
        update_best_solution(inst);

        double t2 = second();

        if(t2 - t1 > inst->time_limit){
            if(VERBOSE>=ERROR){
                print_error("Exceded time limit while computing multi_start_nearest_neighbours, exiting the loop.\n", false);
                break;
            }
        }
    }
}

/**
 * @brief
 * Compute the solution with the variable veighbourhood algorithm, starting with a nearest neighbour solution
 * @param inst the tsp instance
 */
void variable_neighbourhood(instance *inst){
    //multi_start_nearest_neighbours(inst);
    choose_rand_sol(inst);
    two_opt(inst);
    instance temp_inst;

    int iterations_without_improvement = 0;

    while (second() - inst->t_start < inst->time_limit &&
           iterations_without_improvement < MAX_NO_IMPROVEMENT) {
        
        initialize_instance(&temp_inst);
        temp_inst.nnodes = inst->nnodes;

        allocate_instance(&temp_inst);
        copy_instance(inst, &temp_inst);

        temp_inst.solution_cost = inst->best_cost;
        temp_inst.time_limit = inst->time_limit;
        temp_inst.t_start = inst->t_start;
        
        for (int i = 0; i < KICK; i++) {
            three_opt(&temp_inst);
        }

        temp_inst.solution_cost = compute_solution_cost(&temp_inst, temp_inst.solution);

        two_opt(&temp_inst);
        
        if (temp_inst.solution_cost < inst->best_cost) {
            memcpy(inst->solution, temp_inst.solution, (inst->nnodes + 1) * sizeof(int));
            inst->solution_cost = temp_inst.solution_cost;
            check_solution(inst, false);
            update_best_solution(inst);
            iterations_without_improvement = 0;
        } else {
            iterations_without_improvement++;
        }
        
        free_instance(&temp_inst);
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

    inst->solution_cost = compute_solution_cost(inst, inst->solution);
    check_solution(inst, false);
    update_best_solution(inst);

    free(inserted);
}

/**
 * @brief
 * GRASP algorithm + Two Opt local search
 * @param alpha is the hyperparameter that controls the randomness of the solution
 */
void grasp(instance *inst, double alpha) {
    
	
}

/**
 * @brief
 * Compute the solution with the nearest neighbour heuristic algorithm analyzing all possible starting nodes
 * with respect to the time limit
 */
void all_grasp(instance *inst, double alpha) {
    
}