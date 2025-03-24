#define _CRT_SECURE_NO_WARNINGS
#include <utils.h>
#include <parsers.h>

/**
 * @brief
 * Prints an error message
 * @param err the error string
 * @param terminate true if we want the program to terminate its execution
 */
void print_error(const char *err, bool terminate){
    printf("\n\n ERROR: %s \n\n", err); 
    fflush(NULL);
    
    if(terminate){
        exit(EXIT_FAILURE);
    }
}

/**
 * @brief
 * Allocate the memory for the solutions struct
 * @param sol the solutions struct
 */
void allocate_solution_struct(solutions *sol){
    sol->capacity = 16;
    sol->size = 0;
    sol->all_costs = (double*)malloc(sol->capacity * sizeof(double));
}

/**
 * @brief
 * Free the memory allocated for the solutions
 * @param sol the solutions instance
 */
void free_solution_struct(solutions *sol){
    free(sol->all_costs);
}

/**
 * @brief
 * Add a new solution to the solutions struct and eventually realloc dynamically the memory
 * @param sol the solutions struct
 * @param cost the new cost of the best solution to add to the struct
 */
void add_solution(solutions *sol, double cost){
    if(sol->size == sol->capacity){
        sol->capacity *= 2;
        sol->all_costs = (double*)realloc(sol->all_costs, sol->capacity * sizeof(double));
    }

    sol->all_costs[sol->size++] = cost;
}

/**
 * @brief 
 * Plot the solution with gnuplot
 * @param inst the tsp instance
 * @param best If true, plot the best_solution; else plot solution
 * @param wait If true, wait for a key press with getchar(); else skip
 */
void plot_solution(instance *inst, bool best){
    #ifdef _WIN32
		FILE *gnuplotPipe = _popen("gnuplot -persistent", "w");
	#else
		FILE *gnuplotPipe = popen("gnuplot", "w");
	#endif

	int *solution = best ? inst->best_solution : inst->solution;
	double cost = best ? inst->best_cost : inst->solution_cost;

	if(solution == NULL)
		print_error("Solution is not initialized", true);

    fprintf(gnuplotPipe, "set terminal qt title 'TSP Solution'\n");
	fprintf(gnuplotPipe, "set title 'Algorithm: %s, Solution Cost: %.4lf, Time Limit: %.2lf'\n", get_algorithm_name(inst->algorithm), cost, inst->time_limit);
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set grid\n");
	fprintf(gnuplotPipe, "set key outside top\n");

	fprintf(gnuplotPipe, "plot '-' with lines linecolor 'gray' linewidth 2 title 'Edges', '-' with points pointtype 7 pointsize 1.5 linecolor 'blue' title 'Nodes', '-' with points pointtype 7 pointsize 1.5 linecolor 'red' title 'Starting Node'\n");

	for(int i = 0; i < inst->nnodes + 1; i++){
        int idx = inst->solution[i];
        fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
    }
    fprintf(gnuplotPipe, "e\n");

	for(int i = 0; i < inst->nnodes + 1; i++){
        int idx = inst->solution[i];
        fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
    }

    fprintf(gnuplotPipe, "e\n");

	fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[0], inst->ycoord[0]);
	fprintf(gnuplotPipe, "e\n");

	fflush(gnuplotPipe);

    fprintf(gnuplotPipe, "pause mouse close\n");

    #ifdef _WIN32
		_pclose(gnuplotPipe);
	#else
		pclose(gnuplotPipe);
	#endif
}

void plot_solutions(instance *inst){
   #ifdef _WIN32
		FILE *gnuplotPipe = _popen("gnuplot -persistent", "w");
	#else
		FILE *gnuplotPipe = popen("gnuplot", "w");
	#endif

    fprintf(gnuplotPipe, "set terminal pngcairo size 1920,1080 enhanced font 'Verdana,12'\n");
    fprintf(gnuplotPipe, "set output 'history_plot.png'\n");
    fprintf(gnuplotPipe, "set xlabel 'Iteration'\n");
    fprintf(gnuplotPipe, "set ylabel 'Cost'\n");
    fprintf(gnuplotPipe, "set grid\n");
	fprintf(gnuplotPipe, "set key outside top\n");

	fprintf(gnuplotPipe, "plot '-' with lines linecolor 'red' linewidth 2 title 'Best Costs', '-' with lines linecolor 'blue' linewidth 2 title 'Solution Costs'\n");

	for(int i = 0; i < inst->history_best_costs.size; i++){
        fprintf(gnuplotPipe, "%d %lf\n", i, inst->history_best_costs.all_costs[i]);
    }
    fprintf(gnuplotPipe, "e\n");

    for(int i = 0; i < inst->history_costs.size; i++){
        fprintf(gnuplotPipe, "%d %lf\n", i, inst->history_costs.all_costs[i]);
    }
    fprintf(gnuplotPipe, "e\n");

	fflush(gnuplotPipe);
    fprintf(gnuplotPipe, "unset output\n");

    #ifdef _WIN32
		_pclose(gnuplotPipe);
	#else
		pclose(gnuplotPipe);
	#endif
}

/**
 * @brief 
 * Choose which algorithm must be used to solve the problem
 * @param inst the tsp instance
 */
void choose_run_algorithm(instance *inst){
	double t2;

    if(VERBOSE >= INFO)
	    printf("\nMaximum time to solve this problem: %lf seconds\n", inst->time_limit);
    
    inst->t_start = second();

    switch(inst->algorithm){
        case 'N':
            printf("Solving problem with nearest neighbour algorithm\n");
            multi_start_nearest_neighbours(inst);
            break;
        case 'E':
            printf("Solving problem with extra mileage algorithm\n");
            extra_mileage(inst);
            break;
        case 'V':
            if(inst->time_limit == INF_COST){
                printf("Please, insert time limit (seconds): ");
                scanf("%lf", &inst->time_limit);
                getchar();
                printf("Updated time to solve this problem: %lf seconds\n", inst->time_limit);
		    }

            printf("Solving problem with variable neighbourhood algorithm\n");
            variable_neighbourhood(inst);
            break;
        case 'G':
            printf("Solving problem with GRASP\n");
            multi_start_grasp(inst);
            break;
        case 'T':
            if(inst->time_limit == INF_COST){
                printf("Please, insert time limit (seconds): ");
                scanf("%lf", &inst->time_limit);
                getchar();
                printf("Updated time to solve this problem: %lf seconds\n", inst->time_limit);
		    }

            printf("Solving problem with tabu search algorithm\n");
            tabu(inst);
            break;
        default:
            print_error("Algorithm is not available\n", true);
            break;
    }

	t2 = second();
	if (VERBOSE >= INFO)
		printf("\nTSP problem solved in %lf sec.s\n", t2-inst->t_start);
    
	plot_solution(inst, true);
    plot_solutions(inst);
}

/**
 * @brief 
 * Run the same algorithm several time with different time limit on the same file
 * @param inst the tsp instance
 */
void benchmark_algorithm_by_time(instance *inst){
    char algorithm;
    double timeLimits[3];
    double bestCosts[3];
    int *bestSolutions[3] = { NULL, NULL, NULL };

    printf("Choose the algorithm to use for benchmarking:\n");
    printf("N = nearest neighbour\nE = extra mileage\nV = variable neighbourhood\n");
    algorithm = toupper(getchar());
    getchar();

    printf("Insert 3 different time limits (in seconds):\n");
    for(int i = 0; i < 3; i++){
        printf("Time limit %d: ", i + 1);
        scanf("%lf", &timeLimits[i]);
    }
    getchar();

    if(algorithm != 'N' && algorithm != 'E' && algorithm != 'V'){
        printf("Algorithm %c is not available\n", algorithm);
        exit(EXIT_FAILURE);
    }

    for(int i = 0; i < 3; i++){
        printf("\n=== Benchmark with time limit = %lf seconds ===\n", timeLimits[i]);

        inst->time_limit = timeLimits[i];
        inst->t_start = second();

        switch(algorithm){
            case 'N':
                printf("Running multi_start_nearest_neighbours...\n");
                multi_start_nearest_neighbours(inst);
                break;
            case 'E':
                printf("Running extra_mileage...\n");
                extra_mileage(inst);
                break;
            case 'V':
                printf("Running variable_neighbourhood...\n");
                variable_neighbourhood(inst);
                break;
            default:
                printf("Algorithm %c is not available\n", algorithm);
                exit(EXIT_FAILURE);
        }

        double t2 = second();
        bestCosts[i] = inst->best_cost;
        bestSolutions[i] = (int*)calloc(inst->nnodes + 1, sizeof(int));
        memcpy(bestSolutions[i], inst->best_solution, (inst->nnodes + 1) * sizeof(int));
        
    }

    printf("\n=== Benchmark completed. Summary of best costs: ===\n");
    for(int i = 0; i < 3; i++){
        printf("Time limit: %lf, best cost found: %lf\n", timeLimits[i], bestCosts[i]);
    }
    printf("===========================================\n");

    printf("\nDo you want to plot the 3 best solutions? (y/n): ");
    char choice = toupper(getchar());
    getchar();

    if(choice == 'Y'){
        for(int i = 0; i < 3; i++){
            memcpy(inst->best_solution, bestSolutions[i], (inst->nnodes + 1) * sizeof(int));
            inst->time_limit = timeLimits[i];
            inst->best_cost = bestCosts[i];

            printf("\nPlotting solution for time limit = %lf\n", timeLimits[i]);
            plot_solution(inst, true);
        }
    }

    for(int i = 0; i < 3; i++){
        if(bestSolutions[i]){
            free(bestSolutions[i]);
            bestSolutions[i] = NULL;
        }
    }

    printf("\nEnd of benchmark.\n");
}

/**
 * @brief 
 * Function that runs the same algorithm on NUM_FILES different files, and store the best solution cost in a csv file
 */
void benchmark_algorithm_by_params()
{   
    char algorithm;
    double timeLimit;
    //char fileNames[MAX_ROWS - 1][256];
    double bestCosts[MAX_ROWS - 1];

    printf("Choose the algorithm to use for benchmarking:\n");
    printf("N = nearest neighbour\nE = extra mileage\nV = variable neighbourhood\nG = GRASP\nT = Tabu Search\n");
    algorithm = toupper(getchar());
    getchar();

    printf("Enter the time limit (seconds) for all %d runs: ", MAX_ROWS - 1);
    scanf("%lf", &timeLimit);
    while (getchar() != '\n');

    /*for (int i = 0; i < MAX_ROWS - 1; i++) {
        printf("Enter path for file #%d: ", i + 1);
        scanf("%255s", fileNames[i]); 
    }
    while (getchar() != '\n');*/

    for (int i = 0; i < MAX_ROWS - 1; i++) {
        instance inst;

        initialize_instance(&inst);
        //strcpy(inst.input_file, fileNames[i]);
        strcpy(inst.input_file, staticFileNames[i]);
        read_input(&inst);
        compute_all_costs(&inst);
        inst.time_limit = timeLimit;
        inst.t_start = second();

        switch (algorithm) {
            case 'N':
                printf("Running nearest neighbour on %s...\n", staticFileNames[i]);
                multi_start_nearest_neighbours(&inst);
                break;
            case 'E':
                printf("Running extra mileage on %s...\n", staticFileNames[i]);
                extra_mileage(&inst);
                break;
            case 'V':
                printf("Running variable neighbourhood on %s...\n", staticFileNames[i]);
                variable_neighbourhood(&inst);
                break;
            case 'G':
                printf("Running GRASP on %s...\n", staticFileNames[i]);
                multi_start_grasp(&inst);
                break;
            case 'T':
                printf("Running tabu search on %s...\n", staticFileNames[i]);
                tabu(&inst);
                break;
            default:
                printf("Algorithm %c is not available\n", algorithm);
                exit(EXIT_FAILURE);
        }
        bestCosts[i] = inst.best_cost;
        free_instance(&inst);

    }

    char algorithmID[64];
    snprintf(algorithmID, sizeof(algorithmID), "%c_%.0fs", algorithm, timeLimit);
    write_csv(bestCosts, algorithmID);

}