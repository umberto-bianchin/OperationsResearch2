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
	strcpy(inst->input_file, "NULL");

	inst->xcoord = NULL;
	inst->ycoord = NULL;
	inst->best_solution.path = NULL;
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
 * --------------- ATTENTO CHE QUESTO COPIA LA BEST SOLUTION ----------------- DA RIGUARDARE
 * @brief 
 * Copy the parameters of an instance into a new one
 * @param inst the tsp instance to copy
 * @param new_inst the new tsp instance
 */
void copy_instance(instance *inst, instance *new_inst){
	// memcpy(new_inst->solution, inst->best_solution, (inst->nnodes + 1) * sizeof(int));
	memcpy(new_inst->costs, inst->costs, inst->nnodes * inst->nnodes * sizeof(double));
	memcpy(new_inst->xcoord, inst->xcoord, inst->nnodes * sizeof(double));
	memcpy(new_inst->ycoord, inst->ycoord, inst->nnodes * sizeof(double));
	memcpy(new_inst->best_solution.path, inst->best_solution.path, (inst->nnodes + 1) * sizeof(int));
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
	add_solution(&(inst->history_costs), s->cost);
	add_solution(&(inst->history_best_costs), inst->best_solution.cost);
}

/**
 * @brief 
 * Compute the costs of all the edges of the instance
 * @param inst the tsp instance
 */
void compute_all_costs(instance *inst){	
	if(inst->costs == NULL)
		print_error("The costs vector is not initialize\n");

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

	if(fabs(calculated_cost - s->cost) > EPS_COST){
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
	if(s->cost >= inst->best_solution.cost)
		return;	

	check_solution(inst, s);

	if(VERBOSE >= DEBUG)
		printf("Best solution updated: best cost was %f, now is %f\n", inst->best_solution.cost, s->cost);
	
	inst->best_solution.cost = s->cost;

	for(int i = 0; i < inst->nnodes + 1; i++)
		inst->best_solution.path[i] = s->path[i];
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
	
	return round(dis);
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
 * DA AGGIUSTARE, GLI DEVO PASSARE i + 1 e j, ERRORE
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
void two_opt(instance *inst){
	bool improved = true;
	double elapsed_time = second() - inst->t_start;
	solution s;
	copy_solution(&s, &inst->best_solution, inst->nnodes);

	while (improved){	
		improved = false;
		double min_delta = INF_COST;
		int swap_i = -1, swap_j = -1;
		
		for (int i = 1; i < inst->nnodes; i++) {
			for (int j = i + 2; j < inst->nnodes; j++) {

				double current_delta = calculate_delta(i, j, inst, &s);

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

			reverse_segment(swap_i, swap_j, &s);
			compute_solution_cost(inst, &s);
			improved = true;
			elapsed_time = second() - inst->t_start;

			if(elapsed_time > inst->time_limit){
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
 * @param a first node
 * @param b second node
 * @param c third node
 * @param d fourth node
 * @param e fifth node
 * @param f sixth node
 * @param n seventh node
 * @return int best case
 */
int find_best_move(instance *inst, int a, int b, int c, int d, int e, int f, int n){
    double ab = inst->costs[a * inst->nnodes + b], cd = inst->costs[c * inst->nnodes + d], ef = inst->costs[e * inst->nnodes + f];
    double ae = inst->costs[a * inst->nnodes + e], bf = inst->costs[b * inst->nnodes + f];
    double df = inst->costs[d * inst->nnodes + f], ac = inst->costs[a * inst->nnodes + c];
    double be = inst->costs[b * inst->nnodes + e], ad = inst->costs[a * inst->nnodes + d];
	double ec = inst->costs[e * inst->nnodes + c], db = inst->costs[d * inst->nnodes + b];
	double cf = inst->costs[c * inst->nnodes + f], eb = inst->costs[e * inst->nnodes + b];

	double base = ab + cd + ef;
    double gains[4] = { 
        ac + be + df - base, 
        ae + db + cf - base, 
        ad + ec + bf - base, 
        ad + eb + cf - base 
    };

    double maxGain = 0;
	int bestCase = 0;
    for (int i = 0; i < 4; i++) {
        if (gains[i] < 0 && gains[i] < maxGain) {
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
void three_opt(instance *inst){
	double elapsed_time = second() - inst->t_start;
	int i, j, k, temp;
	int nodes = inst->nnodes;
	solution s;
	copy_solution(&s, &inst->best_solution, nodes);

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
	
	int move = find_best_move(inst, s.path[i], s.path[i+1], s.path[j], s.path[j+1], s.path[k], s.path[k+1], nodes);
	
	if(VERBOSE >= DEBUG)
		plot_solution(inst, false);
	
	apply_best_move(inst, i, j, k, move, &s);
}