#define _CRT_SECURE_NO_WARNINGS
#include <tsp_utilities.h>
#include <heuristics.h>

void initialize_instance(instance *inst){
	inst->best_cost = INF_COST;
	inst->solution_cost = INF_COST;
	inst->time_limit = INF_COST;
	inst->seed = 0;
	inst->nnodes = -1;
	inst->t_start = second();
	strcpy(inst->input_file, "NULL");

	inst->xcoord = NULL;
	inst->ycoord = NULL;
	inst->best_solution = NULL;
	inst->solution = NULL;
	inst->costs = NULL;
}

void allocate_instance(instance *inst){
	inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	 
	inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
	inst->best_solution = (int *) calloc(inst->nnodes + 1, sizeof(int));
	inst->solution = (int *) calloc(inst->nnodes + 1, sizeof(int));  
	inst->costs = (double *) calloc(inst->nnodes * inst->nnodes, sizeof(double));	
}

/**
 * @brief
 * Free the memory allocated for the instance
 */
void free_instance(instance *inst){     
	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->costs);
	free(inst->solution);
	free(inst->best_solution);
}

/**
 * @brief 
 * Chooses a valid random solution
 */
void choose_rand_sol(instance *inst){
    srand(inst->seed);
	bool *choosen = calloc(inst->nnodes, sizeof(bool));
	int node;

	for (int i = 0; i < inst->nnodes; i++) {
		do{
			node = rand() % inst->nnodes;
		}while(choosen[node]);

       	inst->solution[i] = node;
		choosen[node] = true;
    }
	free(choosen);

	inst->solution[inst->nnodes] = inst->solution[0];
	inst->solution_cost = compute_solution_cost(inst, inst->solution);
	check_solution(inst, false);
    if(VERBOSE >= DEBUG){
        printf("Choosen value for best_solution: ");
        for (int i = 0; i < inst->nnodes + 1; i++)
			printf("%d ", inst->solution[i]);
		printf("\n");
    }
}

/**
 * @brief 
 * Plot the solution with gnuplot
 * If best = 1 it plots the best solution, otherwise the current solution
 */
void plot_solution(instance *inst, char best){
    #ifdef _WIN32
		FILE *gnuplotPipe = _popen("gnuplot -persistent", "w");
	#else
		FILE *gnuplotPipe = popen("gnuplot ", "w");
	#endif

	int *solution = best ? inst->best_solution : inst->solution;

	if(solution == NULL)
		print_error("Solution is not initialized", true);

	fprintf(gnuplotPipe, "set title 'TSP Solution'\n");
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set grid\n");
	fprintf(gnuplotPipe, "set key outside top\n");
	fprintf(gnuplotPipe, "set term qt title 'TSP Solution'\n");
    fprintf(gnuplotPipe, "plot '-' with linespoints linestyle 1 linewidth 2 pointtype 7 pointsize 1.5 linecolor 'blue' title 'TSP Solution'\n");

	for(int i = 0; i < inst->nnodes + 1; i++)
	{
		int idx = inst->solution[i];
		fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
	}

    fprintf(gnuplotPipe, "e\n");

	fflush(gnuplotPipe);

	getchar();

    #ifdef _WIN32
		_pclose(gnuplotPipe);
	#else
		pclose(gnuplotPipe);
	#endif
}


double compute_solution_cost(instance *inst, int *tour){
    double total_cost = 0.0;
    for (int i = 0; i < inst->nnodes; i++)
        total_cost += inst->costs[tour[i] * inst->nnodes + tour[i + 1]];
    
    return total_cost;
}

/**
 * @brief 
 * Compute the costs of all the edges of the graph
 */
void compute_all_costs(instance *inst){	
	if(inst->costs == NULL)
		print_error("The costs vector is not initialize\n", true);

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
 * Stops the program if the solution is not valid.
 * A solution is valid if:
 * - the first element is equal to the last one
 * - each element is present only once
 * - the solution contains only valid nodes
 * - the cost of the solution is correct
 * 
 * @param best if true checks the best solution, otherwise the current solution
 */
void check_solution(instance *inst, bool best){	
	char error = 0;
	int *solution = best ? inst->best_solution : inst->solution;

	// checks if the first element is equal to the last one
	if(solution[0] != solution[inst->nnodes]) {
		if(VERBOSE >= INFO)
			printf("First element is not equal to the last one.\n");
		
		if(VERBOSE >= DEBUG) 
			printf("First element: %d, last element: %d\n", solution[0], solution[inst->nnodes]);

		free_instance(inst);
		print_error("Solution is not valid.\n", true);
	}

	// checks if each element is present only once
	// checks if the solution contains only valid nodes
	int *count = calloc(inst->nnodes, sizeof(int));
	for(int i = 0; i < inst->nnodes; i++){
		count[solution[i]]++;
		if(count[solution[i]] > 1 && solution[i] != 0){
			if(VERBOSE >= INFO) 
				printf("Element %d is present more than once.\n", solution[i]);

			error = 1;
		}
		if(solution[i] < 0 || solution[i] >= inst->nnodes) {
			if(VERBOSE >= INFO) 
				printf("Element %d is not a valid node.\n", solution[i]);
			
			error = 1;
		}
		if(error){
			free(count);
			free(inst);
			print_error("Solution is not valid.", true);
		}
	}
	free(count);

	// checks if the cost of the solution is correct
	double calculated_cost = 0.0;
	for(int i = 0; i < inst->nnodes; i++)
		calculated_cost += inst->costs[solution[i] * inst->nnodes + solution[i + 1]];

	if(fabs(calculated_cost - (best ? inst->best_cost : inst->solution_cost)) > EPS_COST){
		if(VERBOSE >= INFO) 
			printf("Cost of the solution is not correct.\n");
		
		if(VERBOSE >= DEBUG)
			printf("Calculated cost: %lf, solution cost: %lf\n", calculated_cost, (best ? inst->best_cost : inst->solution_cost));

		free(inst);
		print_error("Solution is not valid.", true);
	}
}

/**
 * @brief
 * Updates the best solution if the current solution is better
 */
void update_best_solution(instance *inst){
	// check if the current solution is worst than the best one
	if(inst->solution_cost >= inst->best_cost)
		return;

	if(VERBOSE >= DEBUG)
		printf("Best solution updated: best cost was %f, now is %f\n", inst->best_cost, inst->solution_cost);
	
	inst->best_cost = inst->solution_cost;

	for(int i = 0; i < inst->nnodes + 1; i++)
		inst->best_solution[i] = inst->solution[i];

	check_solution(inst, true);
}

/**
 * @brief 
 * Calculates the euclidean distance between two points of the graphs
 */
double dist(int i, int j, instance *inst){
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	double dis = sqrt(dx*dx+dy*dy);
	
	if (!inst->integer_costs) return dis;
	return round(dis);
}

/**
 * @brief 
 * Calculates the delta of the cost when two edges are canceled and two are created
 */
double calculate_delta(int i, int j, instance *inst) {
	int node_i = inst->solution[i];
	int node_i1 = inst->solution[i+1];
	int node_j = inst->solution[j];
	int node_j1 = inst->solution[j+1];
    
    double delta = (inst->costs[node_i * inst->nnodes + node_j] + inst->costs[node_i1 * inst->nnodes + node_j1]) 
                   - (inst->costs[node_i * inst->nnodes + node_i1] + inst->costs[node_j * inst->nnodes + node_j1]);
    
    return delta;
}


/**
 * @brief 
 * Reverse a segment of edges
 */
void reverse_segment(int start, int end, instance *inst){
    while (start < end) {
        int temp = inst->solution[start];
        inst->solution[start] = inst->solution[end];
        inst->solution[end] = temp;
        start++;
        end--;
    }
}

/**
 * @brief 
 * Refinement method used to trying to improve the current solution without changing the starting node
 */
void two_opt(instance *inst){
	bool improved = true;
	double elapsed_time = second() - inst->t_start;

	while (improved){	
		improved = false;
		double min_delta = INF_COST;
		int swap_i = -1, swap_j = -1;
		
		for (int i = 1; i < inst->nnodes; i++) {
			for (int j = i + 2; j < inst->nnodes; j++) {

				double current_delta = calculate_delta(i, j, inst);

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

			reverse_segment(swap_i + 1, swap_j, inst);
			inst->solution_cost = compute_solution_cost(inst, inst->solution);
			improved = true;
			elapsed_time = second() - inst->t_start;

			if(elapsed_time > inst->time_limit){
				if(VERBOSE>=INFO)
					print_error("Exceded time limit while computing 2-opt, exiting the loop\n", false);

				improved = false;
				break;
			}
		}
	}
}

double find_best_move(instance *inst, int a, int b, int c, int d, int e, int f, int n, double *totCost){

    double ab = inst->costs[a * inst->nnodes + b], cd = inst->costs[c * inst->nnodes + d], ef = inst->costs[e * inst->nnodes + f];
    double ae = inst->costs[a * inst->nnodes + e], bf = inst->costs[b * inst->nnodes + f], ce = inst->costs[c * inst->nnodes + e];
    double df = inst->costs[d * inst->nnodes + f], ac = inst->costs[a * inst->nnodes + c], bd = inst->costs[b * inst->nnodes + d];
    double be = inst->costs[b * inst->nnodes + e], ad = inst->costs[a * inst->nnodes + d], ec = inst->costs[e * inst->nnodes + c];
    double db = inst->costs[d * inst->nnodes + b], cf = inst->costs[c * inst->nnodes + f], eb = inst->costs[e * inst->nnodes + b];

    double gains[8] = { 
        0, 
        ae + bf - ab - ef, 
        ce + df - cd - ef, 
        ac + bd - ab - cd, 
        ac + be + df - (ab + cd + ef), 
        ae + db + cf - (ab + cd + ef), 
        ad + ec + bf - (ab + cd + ef), 
        ad + eb + cf - (ab + cd + ef) 
    };

    double maxGain = 0, bestCase = 0;
    for (int i = 1; i < 8; i++) {
        if (gains[i] < 0 && gains[i] < maxGain) {
            maxGain = gains[i];
            bestCase = i;
        }
    }

    *totCost += maxGain;
    return bestCase;
}

void apply_best_move(instance *inst, int i, int j, int k, int best_case){
    switch (best_case) {
        case 1:
            reverse_segment(i + 1, j, inst);
            break;
        case 2:
            reverse_segment(j + 1, k, inst);
            break;
        case 3:
            reverse_segment(i + 1, k, inst);
            break;
        case 4:
            reverse_segment(i + 1, j, inst);
            reverse_segment(j + 1, k, inst);
            break;
        case 5:
            reverse_segment(i + 1, j, inst);
            reverse_segment(i + 1, k, inst);
            break;
        case 6:
            reverse_segment(j + 1, k, inst);
            reverse_segment(i + 1, k, inst);
            break;
        case 7:
            reverse_segment(i + 1, j, inst);
            reverse_segment(j + 1, k, inst);
            reverse_segment(i + 1, k, inst);
            break;
        default:
            break;
    }
}



void three_opt(instance *inst){
	bool improved = true;
	double elapsed_time = second() - inst->t_start;

    while (improved) {
        improved = false;
        int best_i = -1, best_j = -1, best_k = -1;
        double best_cost = inst->solution_cost;
        int best_case = 0;

        for (int i = 0; i < inst->nnodes; i++) {
            for (int j = i + 2; j < inst->nnodes; j++) {
                for (int k = j + 2; k < inst->nnodes; k++) {
                    int a = inst->solution[i];
                    int b = inst->solution[i + 1];
                    int c = inst->solution[j];
                    int d = inst->solution[j + 1];
                    int e = inst->solution[k];
                    int f = inst->solution[(k + 1) % inst->nnodes];
                    
                    double new_cost = inst->solution_cost;
                    int move = find_best_move(inst, a, b, c, d, e, f, inst->nnodes, &new_cost);
                    
                    if (new_cost < best_cost) {
                        best_cost = new_cost;
                        best_i = i;
                        best_j = j;
                        best_k = k;
                        best_case = move;
                	}
				}
            }
        }

		if (best_cost < inst->solution_cost) {
            apply_best_move(inst, best_i, best_j, best_k, best_case);

            inst->solution_cost = best_cost;
            improved = 1;
        }

		elapsed_time = second() - inst->t_start;

		if(elapsed_time > inst->time_limit){
			if(VERBOSE>=INFO)
				print_error("Exceded time limit while computing 3-opt, exiting the loop\n", false);
			
			improved = false;
			break;
		}
    }
}


/**
 * @brief 
 * Choose which algorithm must be used to solve the problem
 */
void choose_run_algorithm(instance *inst){
	char algorithm;
	double t1, t2;

	printf("Choose the algorithm to use: N for nearest neighbour, E for extra-mileage, V for variable neighborhood\n");
	algorithm = getchar();
	getchar();

	t1 = second();
	inst->t_start = second();

	if(algorithm == 'N')
	{	
		printf("Solving problem with nearest neighbour algorithm\n");
		multi_start_nearest_neighbours(inst);

	} else if (algorithm == 'E'){
		printf("Solving problem with extra mileage algorithm\n");
		extra_mileage(inst);

	}else if (algorithm == 'V'){
		printf("Solving problem with variable neighbourhood algorithm\n");
		variable_neighbourhood(inst);

	} else{
		printf("Algorithm %c is not available\n", algorithm);
		exit(EXIT_FAILURE);
	}

	t2 = second();
	if (VERBOSE >= INFO)
		printf("\nTSP problem solved in %lf sec.s\n", t2-t1);
}
