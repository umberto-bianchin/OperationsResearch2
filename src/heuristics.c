#include <heuristics.h>

/**
 * @brief
 * Compute the solution with the nearest neighbour heuristic algorithm given a starting node
 */
void nearest_neighbour(instance *inst, int start_node){
    int *visited = (int *) calloc(inst->nnodes, sizeof(int));
    inst->solution[0] = start_node;
    inst->solution[inst->nnodes] = start_node;

    visited[start_node] = 1;
    
    for(int i = 1; i < inst->nnodes; i++){
        int last_selected_node = inst->solution[i-1];
        int nearest_node = -1;
        double min_cost = INF_COST;

        for(int j = 0; j < inst->nnodes; j++){
            if(visited[j] == 0 && inst->costs[last_selected_node * inst->nnodes + j] < min_cost){
                min_cost = inst->costs[last_selected_node * inst->nnodes + j];
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
 */
void multi_start_nearest_neighbours(instance *inst){
    double t1 = second();
    
    for(int i = 0; i < inst->nnodes; i++){
        nearest_neighbour(inst, i);
        update_best_solution(inst);

        double t2 = second();

        if(t2 - t1 > inst->time_limit){
            if(VERBOSE>=INFO){
                print_error("Exceded time limit while computing all_nearest_neighbours, exiting the loop.\n", false);
                break;
            }
        }
    }
}

/**
 * @brief
 * Variable Neighbourhood Search algorithm
 * @param kick is the number of time that the 3-opt algorithm is called for each local optimum solution
 */
void variable_neighbourhood(instance *inst){
    multi_start_nearest_neighbours(inst);
    instance temp_inst;

    double elapsed_time = second() - inst->t_start;

    while(elapsed_time < inst->time_limit){
        initialize_instance(&temp_inst);
        temp_inst.nnodes = inst->nnodes;
        allocate_instance(&temp_inst);

        memcpy(temp_inst.solution, inst->solution, (inst->nnodes + 1) * sizeof(int));
        memcpy(temp_inst.costs, inst->costs, (inst->nnodes + 1) * sizeof(int));
        temp_inst.solution_cost = inst->solution_cost;
        temp_inst.time_limit = inst->time_limit;
        temp_inst.t_start = inst->t_start;

        for(int k = 0; k < KICK; k++){
            three_opt(&temp_inst);
        }

        two_opt(&temp_inst);

        if (temp_inst.solution_cost < inst->solution_cost){
            memcpy(inst->solution, temp_inst.solution, (inst->nnodes + 1) * sizeof(int));
            inst->solution_cost = temp_inst.solution_cost;

            check_solution(inst, false);
            update_best_solution(inst);
        }

        elapsed_time = second() - inst->t_start;
    }

    free_instance(&temp_inst);
}

/**
 * @brief
 * Compute the solution with the extra mileage heuristic algorithm starting from the most distant points
 */
void extra_mileage(instance *inst){
    int i, j, a, b, h, node_a, node_b, node_h, nInserted = 0; 
    double bestDelta, currentDelta, distance, maxDist; 
    int *inserted = (int *) calloc(inst->nnodes, sizeof(int));

    maxDist = -1;
    for(i = 0; i < inst->nnodes; i++){
        for(j = i + 1; j < inst->nnodes; j++){
            distance = inst->costs[i * inst->nnodes + j];
            if(distance > maxDist){
                maxDist = distance;
                node_a = i;
                node_b = j;
            }
        }
    }

    inst->solution[0] = node_a;
    inst->solution[1] = node_b;
    nInserted = 2;
    inserted[node_a] = 1;
    inserted[node_b] = 1;

    while(nInserted < inst->nnodes) {
        bestDelta = INF_COST;
        node_a = -1; node_b = -1; node_h = -1;
        
        for(int pos = 0; pos < nInserted; pos++) {
            a = inst->solution[pos];
            b = inst->solution[pos+1];
            
            for(h = 0; h < inst->nnodes; h++){
                if(inserted[h] == 0) {
                    currentDelta = inst->costs[a * inst->nnodes + h] + inst->costs[h * inst->nnodes + b] - inst->costs[a * inst->nnodes + b];
                    if(currentDelta < bestDelta){
                        bestDelta = currentDelta;
                        node_a = a;
                        node_b = b;
                        node_h = h;
                    }
                }
            }
        }
    

        if(node_h != -1) {
            int pos = 0;
            for(pos = 0; pos < nInserted; pos++){
                if(inst->solution[pos] == node_a && inst->solution[pos+1] == node_b)
                    break;
            }
            
            for(i = nInserted + 1; i > pos+1; i--){
                inst->solution[i] = inst->solution[i - 1];
            }
            
            inst->solution[pos+1] = node_h;
            inserted[node_h] = 1;
            nInserted++;

        } else {
            break;
        }
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