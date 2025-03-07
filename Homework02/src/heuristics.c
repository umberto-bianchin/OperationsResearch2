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
    two_opt(inst);
    check_solution(inst, 0);
}

/**
 * @brief
 * Compute the solution with the nearest neighbour heuristic algorithm analyzing all possible starting nodes
 */
void all_nearest_neighbours(instance *inst){
    inst->best_cost = INF_COST;
    
    double t1 = second();
    
    for(int i = 0; i < inst->nnodes; i++){
        nearest_neighbour(inst, i);
        update_best_solution(inst);

        double t2 = second();

        if(t2 - t1 > inst->time_limit){
            if(VERBOSE>=20){printf("Exceded time limit while computing all_nearest_neighbours, exiting the loop\n");}
            break;
        }
    }
}