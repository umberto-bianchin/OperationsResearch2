#include "tsp_utilities.h"

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }

void free_instance(instance *inst)
{     
	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->costs);
	free(inst->solution);
	free(inst->best_solution);
}

/**
 * @brief 
 * Chooses a random solution
 */
void choose_rand_sol(instance *inst)
{
    srand(inst->seed);
    for (int i = 0; i < inst->nnodes; i++) {
        inst->solution[i] = rand() % inst->nnodes;
    }

	inst->solution[inst->nnodes] = inst->solution[0];

    if(VERBOSE >= 50) 
    {
        printf("Choosen value for best_solution: ");
        for (int i = 0; i < inst->nnodes; i++) printf("%d ", inst->solution[i]);
		printf("\n");
    }
}

/**
 * @brief 
 * Plot the solution with gnuplot
 * If best = 1 it plots the best solution, otherwise the current solution
 */
void plot_solution(instance *inst, char best)
{
    #ifdef _WIN32
		FILE *gnuplotPipe = _popen("gnuplot -persistent", "w");
	#else
		FILE *gnuplotPipe = popen("gnuplot", "w");
	#endif

	int *solution = best ? inst->best_solution : inst->solution;

	if(solution == NULL){printf("Solution is not initialized"); exit(1);}

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

	printf("Press any key to close GnuPlot window...\n");
    getchar();


    #ifdef _WIN32
		_pclose(gnuplotPipe);
	#else
		pclose(gnuplotPipe);
	#endif
}

void calc_solution_cost(instance *inst)
{
	inst->solution_cost = 0.0;
	for(int i = 0; i < inst->nnodes - 1; i++)
		inst->solution_cost += inst->costs[inst->solution[i] * inst->nnodes * inst->solution[i + 1]];
}

/**
 * @brief 
 * Compute the costs of all the edges of the graph
 */
void compute_all_costs(instance *inst)
{	
	if(inst->costs == NULL){printf("The costs vector is not initialize\n"); exit(1);};

	for (int i = 0; i < inst->nnodes; i++) {
        for (int j = 0; j < inst->nnodes; j++) {
            inst->costs[i * inst->nnodes + j] = dist(i, j, inst);
        }
    }

	if(VERBOSE >= 100)	//print the costs matrix
	{
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
 * - each element is present only once
 * - the first element is equal to the last one
 * - the solution contains only valid nodes
 * - the costs are correct
 */
void check_solution(instance *inst, char best)
{	
	int *solution = best ? inst->best_solution : inst->solution;

	// checks if each element is present only once
	int *count = calloc(inst->nnodes, sizeof(int));
	for(int i = 0; i < inst->nnodes; i++)
	{
		count[solution[i]]++;
		if(count[solution[i]] > 1 && solution[i] != 0)
		{
			if(VERBOSE >= 20) 
				printf("Solution is not valid: element %d is present more than once\n", solution[i]);

			free(count);
            exit(EXIT_FAILURE);
		}
	}
	free(count);

	// checks if the first element is equal to the last one
	if(solution[0] != solution[inst->nnodes]) 
	{
		if(VERBOSE >= 20) 
			printf("Solution is not valid: first element is not equal to the last one\n");
		
		if(VERBOSE >= 100) 
			printf("First element: %d, last element: %d\n", solution[0], solution[inst->nnodes]);

        exit(EXIT_FAILURE);
	}

	// checks if the solution contains only valid nodes
	for(int i = 0; i < inst->nnodes; i++)
	{
		if(solution[i] < 0 || solution[i] >= inst->nnodes) 
		{
			if(VERBOSE >= 20) 
				printf("Solution is not valid: element %d is not a valid node\n", solution[i]);
			
			exit(EXIT_FAILURE);
		}
	}

	// checks if the cost of the solution is correct
	double calculated_cost = 0.0;
	for(int i = 0; i < inst->nnodes - 1; i++)
		calculated_cost += inst->costs[solution[i] * inst->nnodes * solution[i + 1]];

	if(fabs(calculated_cost - (best ? inst->best_cost : inst->solution_cost)) > EPS_COST)
	{
		if(VERBOSE >= 20) 
			printf("Solution is not valid: cost of the solution is not correct\n");
		
		if(VERBOSE >= 100)
			printf("Calculated cost: %lf, solution cost: %lf\n", calculated_cost, (best ? inst->best_cost : inst->solution_cost));

		exit(EXIT_FAILURE);
	}
}

/**
 * @brief
 * Updates the best solution if the current solution is better
 */
void update_best_solution(instance *inst)
{
	// check if the current solution is worst than the best one
	if(inst->solution_cost >= inst->best_cost)
		return;

	inst->best_cost = inst->solution_cost;

	for(int i = 0; i < inst->nnodes; i++)
		inst->best_solution[i] = inst->solution[i];
			
	check_solution(inst, 1);
}

/**
 * @brief 
 * Calculates the euclidean distance between two points of the graphs
 */
double dist(int i, int j, instance *inst)
{
	double dx = inst->xcoord[i] - inst->xcoord[j];
	double dy = inst->ycoord[i] - inst->ycoord[j]; 
	double dis = sqrt(dx*dx+dy*dy);

	if ( !inst->integer_costs ) return dis;
	return round(dis);
}

void refine_opt(instance *inst)
{

}