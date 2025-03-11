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
    check_solution(inst, 0);
    two_opt(inst);
    check_solution(inst, 0);
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
    check_solution(inst, 0);
    update_best_solution(inst);

    free(inserted);
}