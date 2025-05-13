#define _CRT_SECURE_NO_WARNINGS
#include <tsp_utilities.h>
#include <heuristics.h>
#include <utils.h>

/**
 * @brief Initialize the TSP instance structure with default parameter values.
 *
 * Must be called before setting node count and allocating any resources.
 *
 * @param inst Pointer to an 'instance' structure to prepare for use.
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

	pthread_mutex_init(&inst->best_mutex, NULL);

	inst->params = (int *) calloc(PARAMS, sizeof(int)); // uses calloc so all entries start at 0

	// Initilalizing all params with best parameters found
	inst->params[KICK] = 7;
	inst->params[K_OPT] = 5;

	inst->params[ALPHA] = 2;
	inst->params[MIN_COSTS] = 3;

	inst->params[MAX_TENURE] = 900;
	inst->params[MIN_TENURE] = 700;
	inst->params[TENURE_STEP] = 50;
	
	inst->params[WARMUP] = 1;
	inst->params[POSTING] = 1;
	inst->params[DEPTH] = 100;
	inst->params[CONCORDE] = 1;

	inst->params[PROBABILITY] = 50;
	inst->params[FIXEDPROB] = 0;

	inst->params[K_LOCAL_BRANCHING] = 20;
}	

/**
 * @brief Allocate memory based on node count stored in the instance.
 *
 * Should follow a call to initialize_instance() and setting inst->nnodes.
 *
 * @param inst Pointer to TSP instance with valid 'nnodes'.
 */
void allocate_instance(instance *inst){
	inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	// zeroed coords 
	inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
	inst->costs = (double *) calloc(inst->nnodes * inst->nnodes, sizeof(double));	
	allocate_solution_struct(&(inst->history_best_costs));
	allocate_solution_struct(&(inst->history_costs));
	allocate_route(&(inst->best_solution), inst->nnodes);
}

/**
 * @brief Allocate and initialize a solution route for 'nnodes'.
 *
 * @param s      Pointer to a solution struct to set up.
 * @param nnodes Number of nodes (tour length is nnodes+1).
 */
void allocate_route(solution *s, int nnodes){
	s->cost = INF_COST;
	s->path = (int *) calloc(nnodes + 1, sizeof(int));	// +1 for return to start
}

/**
 * @brief Copy a solution from src to dest, allocating space as needed.
 *
 * @param dest   Destination solution to receive data.
 * @param src    Source solution to copy from.
 * @param nnodes Number of nodes (length of path minus one).
 */
void copy_solution(solution *dest, solution *src, int nnodes){
	allocate_route(dest, nnodes);
	dest->cost = src->cost;
	memcpy(dest->path, src->path, (nnodes + 1)*sizeof(int));
}

/**
 * @brief Free memory held by a solution's path.
 *
 * @param s Pointer to solution whose path will be freed.
 */
void free_route(solution *s){
	if(s->path != NULL)
		free(s->path);
}


/**
 * @brief Release all dynamic memory in the TSP instance.
 *
 * Ensures no leaks by freeing arrays and solution histories.
 *
 * @param inst Pointer to instance to clean up.
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
	pthread_mutex_destroy(&inst->best_mutex);
}

/**
 * @brief Assign random 2D coordinates to each node in [0, MAX_COORDINATES).
 *
 * Requires inst->seed and inst->nnodes to be set.
 *
 * @param inst Pointer to instance with seed and nnodes.
 */
void set_random_coord(instance *inst){
	srand(inst->seed);
	for (int i = 0; i < inst->nnodes; i++) {
		inst->xcoord[i] = rand() % MAX_COORDINATES;
		inst->ycoord[i] = rand() % MAX_COORDINATES;
	}
}

/**
 * @brief Generate a random valid tour by shuffling node indices.
 *
 * Marks 'best_solution.path' and computes its cost.
 *
 * @param inst Pointer to instance with nnodes and seed set, path allocated.
 */
void choose_rand_sol(instance *inst){
	srand(inst->seed);

	bool *choosen = calloc(inst->nnodes, sizeof(bool));		// track used nodes, zero-init = false
	int node;

	for (int i = 0; i < inst->nnodes; i++) {
		do{
			node = rand() % inst->nnodes;
		}while(choosen[node]);

       	inst->best_solution.path[i] = node;
		choosen[node] = true;
    }

	free(choosen);

	inst->best_solution.path[inst->nnodes] = inst->best_solution.path[0];	// close tour
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
 * @brief Calculate and record the total cost of a tour.
 *
 * Also updates the best-known solution and logs history.
 *
 * @param inst Pointer to instance with cost matrix.
 * @param s    Solution struct containing a path to evaluate.
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
 * @brief Precompute Euclidean distances for all node pairs.
 *
 * Populates inst->costs as a flattened matrix of size nnodes^2.
 *
 * @param inst Pointer to instance with coords and nnodes.
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
 * @brief Validate a tour for completeness and cost consistency.
 *
 * Checks start/end match, unique visits, valid indices, and cost accuracy.
 * Aborts on failure.
 *
 * @param inst Pointer to instance with cost matrix.
 * @param s    Solution to verify.
 */
void check_solution(instance *inst, solution *s){	
	bool error = false;
	int *solution = s->path;

	// Ensure tour is closed
	if(solution[0] != solution[inst->nnodes]) {
		if(VERBOSE >= INFO)
			printf("First element is not equal to the last one.\n");
		
		if(VERBOSE >= DEBUG) 
			printf("First element: %d, last element: %d\n", solution[0], solution[inst->nnodes]);

		free_instance(inst);
		print_error("Solution is not valid.\n");
	}

	// Count visits and check indices
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

	// Verify cost matches sum of edges
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
 * @brief Thread-safe update of the stored best solution if a better one is found.
 *
 * Acquires a mutex before accessing and modifying shared best_solution data.
 *
 * @param inst Pointer to instance holding the current best (with 'best_mutex').
 * @param s    Candidate solution to compare.
 */
void update_best_solution(instance *inst, solution *s){
	pthread_mutex_lock(&inst->best_mutex);

	// check if the current solution is worst than the best one
	if(s->cost >= inst->best_solution.cost - EPS_COST){
		pthread_mutex_unlock(&inst->best_mutex);
		return;	
	}

	check_solution(inst, s);

	if(VERBOSE >= DEBUG)
		printf("Best solution updated: best cost was %f, now is %f\n", inst->best_solution.cost, s->cost);
	
	inst->best_solution.cost = s->cost;

	memcpy(inst->best_solution.path, s->path, (inst->nnodes + 1) * sizeof(int));

	pthread_mutex_unlock(&inst->best_mutex);
}

/**
 * @brief Compute Euclidean or integer-rounded distance between two nodes.
 *
 * @param i    Index of first node.
 * @param j    Index of second node.
 * @param inst Instance with coordinates and integer flag.
 * @return Distance (rounded if inst->integer_costs != 0).
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
 * @brief Compute cost change for removing edges (i,i+1),(j,j+1) and adding (i,j),(i+1,j+1).
 *
 * Used by 2-opt and tabu to identify beneficial swaps.
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
 * @brief Apply 2-opt refinement until no improvement or timeout.
 *
 * @param inst      TSP instance with time budget.
 * @param s         Current solution to refine.
 * @param timelimit Max allowed CPU time since inst->t_start.
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
 * @brief Evaluate 4 possible reconnections for three-opt and pick the best.
 *
 * @param inst TSP instance with precomputed costs.
 * @param s    Current solution tour to analyze.
 * @param i    Index of first edge (i -> i+1).
 * @param j    Index of second edge.
 * @param k    Index of third edge.
 * @return     Case index [0..3] corresponding to minimal cost delta.
 */
int find_best_move(instance *inst, solution *s, int i, int j, int k){
	int n = inst->nnodes;

	// endpoints of three edges
	int a = s->path[i], b = s->path[i + 1];
	int c = s->path[j], d = s->path[j + 1]; 
	int e = s->path[k], f = s->path[k + 1];

	// compute existing costs
    double ab = inst->costs[a*n + b], cd = inst->costs[c*n + d], ef = inst->costs[e*n + f];
	// compute cross-edge costs for 4 cases
    double ae = inst->costs[a*n + e], bf = inst->costs[b*n + f];
    double df = inst->costs[d*n + f], ac = inst->costs[a*n + c];
    double be = inst->costs[b*n + e], ad = inst->costs[a*n + d];
	double ec = inst->costs[e*n + c], db = inst->costs[d*n + b];
	double cf = inst->costs[c*n + f], eb = inst->costs[e*n + b];

	double base = ab + cd + ef;
	// each entry: new connections minus old base cost
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
 * @brief Apply the chosen three-opt reconnection via reverse_segment.
 *
 * @param inst      TSP instance (unused here).
 * @param i         Index of first breakpoint.
 * @param j         Index of second breakpoint.
 * @param k         Index of third breakpoint.
 * @param best_case Case index from find_best_move().
 * @param s         Solution tour to modify.
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
 * @brief Execute a single random three-opt iteration.
 *
 * Chooses non-overlapping breakpoints, finds best move, and applies it.
 *
 * @param inst TSP instance with seed for rand().
 * @param s    Solution tour to refine.
 */
void three_opt(instance *inst, solution *s){
	int i, j, k, temp;
	int nodes = inst->nnodes;

	// pick three distinct, non-consecutive indices
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
 * @brief Verify nodes selected for k-opt are non-consecutive and unique.
 *
 * @param nodes Array of indices.
 * @param n     Number of nodes in array.
 * @return      true if valid, false otherwise.
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
 * @brief Perform a random k-opt move by reversing k randomly chosen segments.
 *
 * @param inst TSP instance with valid nnodes.
 * @param s    Solution tour to modify.
 * @param k    Number of breakpoints to choose (>=2).
 */
void random_k_opt(instance *inst, solution *s, int k){
	if(k > inst->nnodes || k < 2)
		print_error("Invalid k in k-opt");
	
	int *nodes = calloc(k, sizeof(int));
	int temp;
	
	// select k valid, non-consecutive nodes
	do{
		for(int i = 0; i < k; i++){
			nodes[i] = rand() % inst->nnodes;			
		}
	}while(check_valid_kopt_nodes(nodes, k));

	// sort nodes for consistent segment reversal order
	for(int i = 0; i < k - 1; i++){
		for(int j = i + 1; j < k; j++){
			if(nodes[j] < nodes[i]){
				temp = nodes[i];
				nodes[i] = nodes[j];
				nodes[j] = temp;
			}
		}
	}

	// apply k random 2-opt-like swaps among chosen points
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
 * @brief Apply a random 5-opt move taken by 2 5-opt moves by combining several 2-opt segments.
 *
 * @param inst TSP instance with nnodes >=5.
 * @param s    Solution tour to modify.
 */
void five_opt(instance *inst, solution *s){
	// Each element of the array rapresent an edge, nodes[0] is the edge that starts from s->path[edge[0]] and goes to s->path[edge[0] + 1]
	int *edges = calloc(5, sizeof(int));
	int nnodes = inst->nnodes;
	int temp;

	// selects 5 non consecutive nodes
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

	// randomly pick one of two 5-opt patterns
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