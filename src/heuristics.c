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

    calc_solution_cost(inst);
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
 * Variable Neighborhood Search algorithm
 * @param kick is the number of time that the 3-opt algorithm is called for each local optimum solution
 */
void vns(instance *inst, int kick){
    // choose a random solution
    choose_rand_sol(inst);
    check_solution(inst, false);

    while(second() - inst->t_start < inst->time_limit){
        two_opt(inst);
        check_solution(inst, false);
        update_best_solution(inst);

        for(int k = 0; k < kick; k++){
            int *nodes = (int *) calloc(3, sizeof(int));
            do{
                nodes[0] = rand() % inst->nnodes;
                nodes[1] = rand() % inst->nnodes;
                nodes[2] = rand() % inst->nnodes;
            }while(nodes[0] == nodes[1] || nodes[0] == nodes[2] || nodes[1] == nodes[2]);
            
            
        }
    }
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

    calc_solution_cost(inst);
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