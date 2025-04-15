#define _CRT_SECURE_NO_WARNINGS
#include <tsp_utilities.h>
#include <heuristics.h>
#include <utils.h>

/**
 * @brief 
 * Initialize the instance
 * @param inst the tsp instance to initialize
 */
void initialize_instance(instance *inst){
	inst->best_solution.cost = INF_COST;
	inst->time_limit = INF_COST;
	inst->seed = 0;
	inst->nnodes = -1;
	inst->t_start = second();
	inst->algorithm = ' ';
	inst->integer_costs = 0;
	inst->ncols = 0;
	strcpy(inst->input_file, "NULL");

	inst->xcoord = NULL;
	inst->ycoord = NULL;
	inst->best_solution.path = NULL;
	
	inst->params = (int *) calloc(PARAMS, sizeof(int));

	// Initilalizing all params with best parameters found
	inst->params[KICK] = 5;
	inst->params[K_OPT] = 3;

	inst->params[ALPHA] = 20;
	inst->params[MIN_COSTS] = 3;

	inst->params[MAX_TENURE] = 500;
	inst->params[MIN_TENURE] = 200;
	inst->params[TENURE_STEP] = 50;
	
	inst->params[WARMUP] = 0;

}

/**
 * @brief
 * Allocate the memory for the instance's members
 * @param inst the tsp instance
 */
void allocate_instance(instance *inst){
	inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	 
	inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
	inst->costs = (double *) calloc(inst->nnodes * inst->nnodes, sizeof(double));	
	allocate_solution_struct(&(inst->history_best_costs));
	allocate_solution_struct(&(inst->history_costs));
	allocate_route(&(inst->best_solution), inst->nnodes);
}

void allocate_route(solution *s, int nnodes){
	s->cost = INF_COST;
	s->path = (int *) calloc(nnodes + 1, sizeof(int));
}

void copy_solution(solution *dest, solution *src, int nnodes){
	allocate_route(dest, nnodes);
	dest->cost = src->cost;
	memcpy(dest->path, src->path, (nnodes + 1)*sizeof(int));
}

void free_route(solution *s){
	if(s->path != NULL)
		free(s->path);
}


/**
 * @brief
 * Free the memory allocated for the instance
 * @param inst the tsp instance
 */
void free_instance(instance *inst){ 
	if(inst->xcoord != NULL)
		free(inst->xcoord);
	if(inst->ycoord != NULL)
		free(inst->ycoord);
	if(inst->costs != NULL)
		free(inst->costs);

	free_solution_struct(&(inst->history_best_costs));
	free_solution_struct(&(inst->history_costs));
	free_route(&(inst->best_solution));
}

/**
 * @brief
 * Set random coordinates for the nodes of the instance
 * @param inst the tsp instance
 */
void set_random_coord(instance *inst){
	srand(inst->seed);
	for (int i = 0; i < inst->nnodes; i++) {
		inst->xcoord[i] = rand() % MAX_COORDINATES;
		inst->ycoord[i] = rand() % MAX_COORDINATES;
	}
}

/**
 * @brief 
 * Chooses a valid random solution
 * @param inst the tsp instance
 */
void choose_rand_sol(instance *inst){
	srand(inst->seed);

	bool *choosen = calloc(inst->nnodes, sizeof(bool));
	int node;

	for (int i = 0; i < inst->nnodes; i++) {
		do{
			node = rand() % inst->nnodes;
		}while(choosen[node]);

       	inst->best_solution.path[i] = node;
		choosen[node] = true;
    }

	free(choosen);

	inst->best_solution.path[inst->nnodes] = inst->best_solution.path[0];
	compute_solution_cost(inst, &(inst->best_solution));
	check_solution(inst, &(inst->best_solution));

    if(VERBOSE >= DEBUG){
        printf("Choosen value for best_solution: ");
        for (int i = 0; i < inst->nnodes + 1; i++)
			printf("%d ", inst->best_solution.path[i]);
		printf("\n");
    }
}

/**
 * @brief 
 * Compute the cost of the solution
 * @param tour the solution used to compute the cost
 */
void compute_solution_cost(instance *inst, solution *s) {
	double total_cost = 0.0;
	for (int i = 0; i < inst->nnodes; i++)
		total_cost += inst->costs[s->path[i] * inst->nnodes + s->path[i + 1]];
	
	s->cost = total_cost;
	update_best_solution(inst, s);
	add_solution(&(inst->history_costs), s->cost, -1);
	add_solution(&(inst->history_best_costs), inst->best_solution.cost, -1);
}

/**
 * @brief 
 * Compute the costs of all the edges of the instance
 * @param inst the tsp instance
 */
void compute_all_costs(instance *inst){	
	if(inst->costs == NULL)
		print_error("The costs vector is not initialize\n");

	if(inst->nnodes < 3){
		print_error("The istance must have more than 3 nodes.\n");
	}

	for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            inst->costs[i * inst->nnodes + j] = dist(i, j, inst);
        }
    }

	if(VERBOSE >= DEBUG){ //print the costs matrix
		for (int i = 0; i < inst->nnodes; i++) {
			for (int j = 0; j < inst->nnodes; j++) {
				printf("%8.2f ", inst->costs[i * inst->nnodes + j]);
			}
			printf("\n");
		}
	}
}

/**
 * @brief 
 * Check if the solution of the tsp is valid. If not, it stops the program
 * A solution is valid if:
 * - the first element is equal to the last one
 * - each element is present only once
 * - the solution contains only valid nodes
 * - the cost of the solution is correct
 * @param inst the tsp instance
 */
void check_solution(instance *inst, solution *s){	
	bool error = false;
	int *solution = s->path;

	// checks if the first element is equal to the last one
	if(solution[0] != solution[inst->nnodes]) {
		if(VERBOSE >= INFO)
			printf("First element is not equal to the last one.\n");
		
		if(VERBOSE >= DEBUG) 
			printf("First element: %d, last element: %d\n", solution[0], solution[inst->nnodes]);

		free_instance(inst);
		print_error("Solution is not valid.\n");
	}

	// checks if each element is present only once
	// checks if the solution contains only valid nodes
	int *count = calloc(inst->nnodes, sizeof(int));
	for(int i = 0; i < inst->nnodes; i++){
		count[solution[i]]++;
		if(count[solution[i]] > 1 && solution[i] != 0){
			if(VERBOSE >= INFO) 
				printf("Element %d is present more than once.\n", solution[i]);

			error = true;
		}
		if(solution[i] < 0 || solution[i] >= inst->nnodes) {
			if(VERBOSE >= INFO) 
				printf("Element %d is not a valid node.\n", solution[i]);
			
			error = true;
		}
		if(error){
			free(count);
			free(inst);
			print_error("Solution is not valid.");
		}
	}
	free(count);

	// checks if the cost of the solution is correct
	double calculated_cost = 0.0;
	for(int i = 0; i < inst->nnodes; i++)
		calculated_cost += inst->costs[solution[i] * inst->nnodes + solution[i + 1]];

	if(fabs(calculated_cost - s->cost) > EPS_ERR){
		if(VERBOSE >= INFO) 
			printf("Cost of the solution is not correct.\n");
		
		if(VERBOSE >= DEBUG)
			printf("Calculated cost: %lf, solution cost: %lf\n", calculated_cost, s->cost);

		free(inst);
		print_error("Solution is not valid.");
	}
}

/**
 * @brief
 * Updates the best solution if the current solution is better
 * @param inst the tsp instance
 */
void update_best_solution(instance *inst, solution *s){
	// check if the current solution is worst than the best one
	if(s->cost >= inst->best_solution.cost - EPS_COST)
		return;	

	check_solution(inst, s);

	if(VERBOSE >= DEBUG)
		printf("Best solution updated: best cost was %f, now is %f\n", inst->best_solution.cost, s->cost);
	
	inst->best_solution.cost = s->cost;

	memcpy(inst->best_solution.path, s->path, (inst->nnodes + 1) * sizeof(int));
}

/**
 * @brief 
 * Calculates the euclidean distance between two points of the graphs
 * @param i the index of the first node
 * @param j the index of the second node
 * @param inst the tsp instance
 */
double dist(int i, int j, instance *inst){
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	double dis = sqrt(dx*dx+dy*dy);
	
	if(inst->integer_costs)
		return round(dis);
		
	return dis;
}

/**
 * @brief 
 * Calculates the delta of the cost when two edges are canceled and two are created. Used for the two-opt method
 * @param i the index of the first node
 * @param j the index of the second node
 * @param inst the tsp instance
 * @return double delta
 */
double calculate_delta(int i, int j, instance *inst, solution *s){
	int node_i = s->path[i];
	int node_i1 = s->path[i+1];
	int node_j = s->path[j];
	int node_j1 = s->path[j+1];
    
    double delta = (inst->costs[node_i * inst->nnodes + node_j] + inst->costs[node_i1 * inst->nnodes + node_j1]) 
                   - (inst->costs[node_i * inst->nnodes + node_i1] + inst->costs[node_j * inst->nnodes + node_j1]);
    
    return delta;
}


/**
 * @brief 
 * Reverse a segment of edges, swapping tho nodes
 * @param start the index of the first node (the swap actually starts from i+1)
 * @param end the index of the second node
 * @param inst the tsp instance
 */
void reverse_segment(int start, int end, solution *s){
	int swap_i = start + 1;
    while (swap_i < end) {
        int temp = s->path[swap_i];
        s->path[swap_i] = s->path[end];
        s->path[end] = temp;
        swap_i++;
        end--;
    }
}

/**
 * @brief 
 * Refinement method used to trying to improve the current solution
 * @param inst the tsp instance
 */
void two_opt(instance *inst, solution *s, double timelimit){
	bool improved = true;
	double elapsed_time = second() - inst->t_start;

	while (improved){	
		improved = false;
		double min_delta = INF_COST;
		int swap_i = -1, swap_j = -1;
		
		for (int i = 1; i < inst->nnodes; i++) {
			for (int j = i + 2; j < inst->nnodes; j++) {

				double current_delta = calculate_delta(i, j, inst, s);

				if(current_delta < min_delta){
					min_delta = current_delta;
					swap_i = i;
					swap_j = j;
				}
			}
		}

		if(min_delta < 0){
			if(VERBOSE >= DEBUG)
				printf("Swapping node %d with node %d\n", swap_i, swap_j);

			reverse_segment(swap_i, swap_j, s);
			compute_solution_cost(inst, s);
			improved = true;
			elapsed_time = second() - inst->t_start;

			if(elapsed_time > timelimit){
				if(VERBOSE>=ERROR)
					printf("Exceded time limit while computing 2-opt, exiting the loop\n");

				improved = false;
				break;
			}
		}
	}
}

/**
 * @brief 
 * Finds the best move for the three-opt method
 * @param inst the tsp instance
 * @param s is the solution analized
 * @param i first edge (solution path[i] -> path[i + 1])
 * @param j second edge
 * @param k third edge
 * @param maxGain return value, is the delta applied 
 * @return int best case
 */
int find_best_move(instance *inst, solution *s, int i, int j, int k){
	int n = inst->nnodes;

	int a = s->path[i], b = s->path[i + 1];
	int c = s->path[j], d = s->path[j + 1]; 
	int e = s->path[k], f = s->path[k + 1];

    double ab = inst->costs[a*n + b], cd = inst->costs[c*n + d], ef = inst->costs[e*n + f];
    double ae = inst->costs[a*n + e], bf = inst->costs[b*n + f];
    double df = inst->costs[d*n + f], ac = inst->costs[a*n + c];
    double be = inst->costs[b*n + e], ad = inst->costs[a*n + d];
	double ec = inst->costs[e*n + c], db = inst->costs[d*n + b];
	double cf = inst->costs[c*n + f], eb = inst->costs[e*n + b];

	double base = ab + cd + ef;
    double gains[4] = { 
        ac + be + df - base, 
        ae + db + cf - base, 
        ad + ec + bf - base, 
        ad + eb + cf - base 
    };

    double maxGain = INF_COST;
	int bestCase = -1;
    for (int i = 0; i < 4; i++) {
        if (gains[i] < maxGain) {
            maxGain = gains[i];
            bestCase = i;
        }
    }
	
    return bestCase;
}

/**
 * @brief 
 * Apply the best move for the three-opt method, swapping the corresponding nodes
 * @param inst the tsp instance
 * @param i index of the first node
 * @param j index of the second node
 * @param k index of the third node
 * @param best_case the best move choosend by the find_best_move method
 */
void apply_best_move(instance *inst, int i, int j, int k, int best_case, solution *s){
    switch (best_case) {
        case 0:
            reverse_segment(i, j, s);
            reverse_segment(j, k, s);
            break;
        case 1:
            reverse_segment(i, j, s);
            reverse_segment(i, k, s);
            break;
        case 2:
            reverse_segment(j, k, s);
            reverse_segment(i, k, s);
            break;
        case 3:
            reverse_segment(i, j, s);
            reverse_segment(j, k, s);
            reverse_segment(i, k, s);
            break;
        default:
            break;
    }
}

/**
 * @brief 
 * Apply a three-opt move to the instance
 * @param inst the tsp instance
 */
void three_opt(instance *inst, solution *s){
	int i, j, k, temp;
	int nodes = inst->nnodes;

	do {
		i = rand() % nodes;
		j = rand() % nodes;
		k = rand() % nodes;
	} while (abs(i - j) <= 1 || abs(i - k) <= 1 || abs(j - k) <= 1);

	if (i > j){
		temp = i;
		i = j;
		j = temp;
	}
	if (i > k){
		temp = i;
		i = k;
		k = temp;
	}
	if (j > k){
		temp = j;
		j = k;
		k = temp;
	}
	
	int move = find_best_move(inst, s, i, j, k);
	
	if(VERBOSE >= DEBUG)
		plot_solution(inst, s);
	
	apply_best_move(inst, i, j, k, move, s);
}


/**
 * @brief
 * Checks if the selected nodes for the k-opt are valid nodes
 * @return true if all the nodes are different and not consecutive, otherwise false
 */
bool check_valid_kopt_nodes(int *nodes, int n){
	for(int i = 0; i < n - 1; i++){
		for(int j = i + 1; j < n; j++){
			if(abs(nodes[i] - nodes[j]) <= 1)
				return false;
		}
	}

	return true;
}

/**
 * @brief 
 * Apply a random k-opt move to the instance
 * @param inst the tsp instance
 * @param s solution modified
 * @param k number of swap
 */
void random_k_opt(instance *inst, solution *s, int k){
	if(k > inst->nnodes || k < 2)
		print_error("Invalid k in k-opt");
	
	int *nodes = calloc(k, sizeof(int));
	int temp;
	
	// selects k non consecutive nodes
	do{
		for(int i = 0; i < k; i++){
			nodes[i] = rand() % inst->nnodes;			
		}
	}while(check_valid_kopt_nodes(nodes, k));

	for(int i = 0; i < k - 1; i++){
		for(int j = i + 1; j < k; j++){
			if(nodes[j] < nodes[i]){
				temp = nodes[i];
				nodes[i] = nodes[j];
				nodes[j] = temp;
			}
		}
	}

	int idx1, idx2;
    for (int i = 0; i < k; i++) {
		do{
			idx1 = rand() % k;
			idx2 = rand() % k;
		}while(idx1 == idx2);

        reverse_segment(nodes[idx1], nodes[idx2], s);
    }
		
	if(VERBOSE >= DEBUG)
		plot_solution(inst, s);

	free(nodes);
}

/**
 * @brief 
 * Apply a 5-opt move to the instance
 * @param inst the tsp instance
 * @param s solution modified
 */
void five_opt(instance *inst, solution *s){
	// Each element of the array rapresent an edge, nodes[0] is the edge that starts from s->path[edge[0]] and goes to s->path[edge[0] + 1]
	int *edges = calloc(5, sizeof(int));
	int nnodes = inst->nnodes;
	int temp;

	// selects k non consecutive nodes
	do{
		for(int i = 0; i < 5; i++){
			edges[i] = rand() % nnodes;			
		}
	}while(check_valid_kopt_nodes(edges, 5));

	for(int i = 0; i < 5 - 1; i++){
		for(int j = i + 1; j < 5; j++){
			if(edges[j] < edges[i]){
				temp = edges[i];
				edges[i] = edges[j];
				edges[j] = temp;
			}
		}
	}
	
	int A = edges[0];
	int B = edges[1];
	int C = edges[2];
	int D = edges[3];
	int E = edges[4];

	// Select a random valid 5 opt move 
	if(rand()%2 == 0){
		reverse_segment(A, B, s);
		reverse_segment(C, D, s);
		reverse_segment(B, E, s);
		reverse_segment(A, C, s);
		reverse_segment(D, E, s);		
	}
	else{		
		reverse_segment(A, C, s);
		reverse_segment(B, D, s);
		reverse_segment(C, E, s);
		reverse_segment(D, A, s);
		reverse_segment(E, B, s);  
	}
			
	if(VERBOSE >= DEBUG)
		plot_solution(inst, s);

	free(edges);
}