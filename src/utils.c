#define _CRT_SECURE_NO_WARNINGS
#include <utils.h>
#include <parsers.h>
#include <cplex_utilities.h>

/**
 * @brief
 * Prints an error message
 * @param err the error string
 */
void print_error(const char *err){
    printf("\n\n ERROR: %s \n\n", err); 
    fflush(NULL);
    
    exit(EXIT_FAILURE);
    
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
    sol->iteration_times = (double*)malloc(sol->capacity * sizeof(double));
}

/**
 * @brief
 * Free the memory allocated for the solutions
 * @param sol the solutions instance
 */
void free_solution_struct(solutions *sol){
    if(sol->all_costs != NULL)
        free(sol->all_costs);
    if(sol->iteration_times != NULL)
        free(sol->iteration_times);
}

/**
 * @brief
 * Add a new solution to the solutions struct and eventually realloc dynamically the memory
 * @param sol the solutions struct
 * @param cost the new cost of the best solution to add to the struct
 */
void add_solution(solutions *sol, double cost, double time){
    if(sol->size == sol->capacity){
        sol->capacity *= 2;
        sol->all_costs = (double*)realloc(sol->all_costs, sol->capacity * sizeof(double));
        sol->iteration_times = (double*)realloc(sol->iteration_times, sol->capacity * sizeof(double));
    }

    sol->all_costs[sol->size++] = cost;
    if(time != -1)
        sol->iteration_times[sol->size] = time;
}

/**
 * @brief 
 * Plot the solution with gnuplot
 * @param inst the tsp instance
 * @param best If true, plot the best_solution; else plot solution
 * @param wait If true, wait for a key press with getchar(); else skip
 */
void plot_solution(instance *inst, solution *s){
    #ifdef _WIN32
		FILE *gnuplotPipe = _popen("gnuplot -persistent", "w");
	#else
		FILE *gnuplotPipe = popen("gnuplot", "w");
	#endif

	if(s->path == NULL)
		print_error("Solution is not initialized");

    char runID[64];
    char nameComplete[128];
    char timestamp[20];
    setAlgorithmId(inst, runID);

    time_t t = time(NULL);
    struct tm* tm_info = localtime(&t);

    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d_%H-%M-%S", tm_info);

    snprintf(nameComplete, sizeof(nameComplete), "history/solution%s_%s.png", runID, timestamp);

    fprintf(gnuplotPipe, "set terminal pngcairo size 1920,1080 enhanced font 'Verdana,12'\n");
    fprintf(gnuplotPipe, "set output '%s'\n", nameComplete);
	fprintf(gnuplotPipe, "set title 'Algorithm: %s, Solution Cost: %.4lf, Time Limit: %.2lf'\n", print_algorithm(inst->algorithm), s->cost, inst->time_limit);
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set grid\n");
	fprintf(gnuplotPipe, "set key outside top\n");

	fprintf(gnuplotPipe, "plot '-' with lines linecolor 'gray' linewidth 2 title 'Edges', '-' with points pointtype 7 pointsize 1.5 linecolor 'blue' title 'Nodes', '-' with points pointtype 7 pointsize 1.5 linecolor 'red' title 'Starting Node'\n");

	for(int i = 0; i < inst->nnodes + 1; i++){
        int idx = s->path[i];
        fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
    }
    fprintf(gnuplotPipe, "e\n");

	for(int i = 0; i < inst->nnodes + 1; i++){
        int idx = s->path[i];
        fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
    }

    fprintf(gnuplotPipe, "e\n");

	fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[0], inst->ycoord[0]);
	fprintf(gnuplotPipe, "e\n");

	fflush(gnuplotPipe);
    fprintf(gnuplotPipe, "unset output\n");

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

    char runID[64];
    char nameComplete[128];
    char timestamp[20];
    setAlgorithmId(inst, runID);

    time_t t = time(NULL);
    struct tm* tm_info = localtime(&t);

    strftime(timestamp, sizeof(timestamp), "%Y-%m-%d_%H-%M-%S", tm_info);

    snprintf(nameComplete, sizeof(nameComplete), "history/history%s_%s.png", runID, timestamp);

    fprintf(gnuplotPipe, "set terminal pngcairo size 1920,1080 enhanced font 'Verdana,12'\n");
    fprintf(gnuplotPipe, "set output '%s'\n", nameComplete);
    fprintf(gnuplotPipe, "set title 'Algorithm: %s, Solution Cost: %.4lf, Time Limit: %.2lf'\n", print_algorithm(inst->algorithm), inst->best_solution.cost, inst->time_limit);
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

void plot_cplex_solutions(instance *inst){
    #ifdef _WIN32
         FILE *gnuplotPipe = _popen("gnuplot -persistent", "w");
     #else
         FILE *gnuplotPipe = popen("gnuplot", "w");
     #endif
 
     char runID[64];
     char nameComplete[128];
     char timestamp[20];
     setAlgorithmId(inst, runID);
 
     time_t t = time(NULL);
     struct tm* tm_info = localtime(&t);
 
     strftime(timestamp, sizeof(timestamp), "%Y-%m-%d_%H-%M-%S", tm_info);
 
     snprintf(nameComplete, sizeof(nameComplete), "history/history%s_%s.png", runID, timestamp);
 
     fprintf(gnuplotPipe, "set terminal pngcairo size 1920,1080 enhanced font 'Verdana,12'\n");
     fprintf(gnuplotPipe, "set output '%s'\n", nameComplete);
     fprintf(gnuplotPipe, "set title 'Algorithm: %s, Solution Cost: %.4lf, Time Limit: %.2lf'\n", print_algorithm(inst->algorithm), inst->best_solution.cost, inst->time_limit);
     fprintf(gnuplotPipe, "set xlabel 'Time'\n");
     fprintf(gnuplotPipe, "set ylabel 'Cost'\n");
     fprintf(gnuplotPipe, "set grid\n");
     fprintf(gnuplotPipe, "set key outside top\n");
 
     fprintf(gnuplotPipe, "plot '-' with linespoints linecolor 'blue' linewidth 2 pointtype 7 pointsize 1.5 title 'Solution Cost'\n");
 
     for(int i = 0; i < inst->history_best_costs.size; i++){
         fprintf(gnuplotPipe, "%lf %lf\n", inst->history_best_costs.iteration_times[i], inst->history_best_costs.all_costs[i]);
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
	    printf("\nMaximum time to solve this problem: %.3lf seconds\n", inst->time_limit);
    
    inst->t_start = second();

    switch(inst->algorithm){
        case 'N':
            printf("Solving problem with nearest neighbour algorithm\n");
            multi_start_nearest_neighbours(inst, inst->time_limit);
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
            nearest_neighbour(inst, rand() % inst->nnodes);     //need to initialize and optimize the first solution
            solution s;
            copy_solution(&s, &inst->best_solution, inst->nnodes);
            two_opt(inst, &s, inst->time_limit);

            variable_neighbourhood(inst, inst->time_limit);
            break;
        case 'G':
            printf("Solving problem with GRASP\n");
            multi_start_grasp(inst, inst->time_limit);
            break;
        case 'T':
            if(inst->time_limit == INF_COST){
                printf("Please, insert time limit (seconds): ");
                scanf("%lf", &inst->time_limit);
                getchar();
                printf("Updated time to solve this problem: %lf seconds\n", inst->time_limit);
		    }

            printf("Solving problem with tabu search algorithm\n");
            nearest_neighbour(inst, rand() % inst->nnodes);
            tabu(inst, inst->time_limit);
            break;
        case 'B':
            if(inst->time_limit == INF_COST){
                printf("Please, insert time limit (seconds): ");
                scanf("%lf", &inst->time_limit);
                getchar();
                printf("Updated time to solve this problem: %lf seconds\n", inst->time_limit);
            }
            printf("Solving problem with benders\n");
            TSPopt(inst);
            break;
        default:
            print_error("Algorithm is not available\n");
            break;
    }

	t2 = second();
	if (VERBOSE >= INFO)
		printf("\nTSP problem solved in %lf sec.s\n", t2-inst->t_start);
    
	plot_solution(inst, &(inst->best_solution));
    if(inst->algorithm != 'B'){
        plot_solutions(inst);
    } else{
        plot_cplex_solutions(inst);
    }
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
                multi_start_nearest_neighbours(inst, inst->time_limit);
                break;
            case 'E':
                printf("Running extra_mileage...\n");
                extra_mileage(inst);
                break;
            case 'V':
                printf("Running variable_neighbourhood...\n");
                nearest_neighbour(inst, rand() % inst->nnodes);     //need to initialize and optimize the first solution
                solution s;
                copy_solution(&s, &inst->best_solution, inst->nnodes);
                two_opt(inst, &s, inst->time_limit);
                
                variable_neighbourhood(inst, inst->time_limit);
                break;
            default:
                printf("Algorithm %c is not available\n", algorithm);
                exit(EXIT_FAILURE);
        }

        double t2 = second();
        bestCosts[i] = inst->best_solution.cost;
        bestSolutions[i] = (int*)calloc(inst->nnodes + 1, sizeof(int));
        memcpy(bestSolutions[i], inst->best_solution.path, (inst->nnodes + 1) * sizeof(int));
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
            memcpy(inst->best_solution.path, bestSolutions[i], (inst->nnodes + 1) * sizeof(int));
            inst->time_limit = timeLimits[i];
            inst->best_solution.cost = bestCosts[i];

            printf("\nPlotting solution for time limit = %lf\n", timeLimits[i]);
            plot_solution(inst, &(inst->best_solution));
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
void benchmark_algorithm_by_params(instance *inst)
{   
    double bestCosts[MAX_ROWS - 1];

    for (int i = 0; i < MAX_ROWS - 1; i++) {
        inst->seed = i * 1000; //multiply by 100 to obtain more variations
        inst->best_solution.cost = INF_COST;

        set_random_coord(inst);
        compute_all_costs(inst);
        inst->t_start = second();

        switch (inst->algorithm) {
            case 'N':
                printf("Running nearest neighbour on random coordinates with seed %d...\n", inst->seed);
                multi_start_nearest_neighbours(inst, inst->time_limit);
                break;
            case 'E':
                printf("Running extra mileage on random coordinates with seed %d...\n", inst->seed);
                extra_mileage(inst);
                break;
            case 'V':
                printf("Running variable neighbourhood on random coordinates with seed %d...\n", inst->seed);
                nearest_neighbour(inst, rand() % inst->nnodes);     //need to initialize and optimize the first solution
                solution s;
                copy_solution(&s, &inst->best_solution, inst->nnodes);
                two_opt(inst, &s, inst->time_limit);

                variable_neighbourhood(inst, inst->time_limit);
                break;
            case 'G':
                printf("Running GRASP on random coordinates with seed %d...\n", inst->seed);
                multi_start_grasp(inst, inst->time_limit);
                break;
            case 'T':
                printf("Running tabu search on random coordinates with seed %d...\n", inst->seed);
                nearest_neighbour(inst, rand() % inst->nnodes);     //needed to initialize the first solution
                tabu(inst, inst->time_limit);
                break;
            default:
                printf("Algorithm %c is not available\n", inst->algorithm);
                exit(EXIT_FAILURE);
        }
        bestCosts[i] = inst->best_solution.cost;
    }

    char algorithmID[64];
    
    setAlgorithmId(inst, algorithmID);

    write_csv(bestCosts, algorithmID, inst->algorithm);
}

/**
 * @brief
 * Check if the algorithm is valid, stop the program if is not
 * @param alg struct containing the algorithms
 * @param algorithm the algorithm to check
 */
void check_valid_algorithm(char algorithm){
    bool valid = false;

    for(int i = 0; i < ALGORITHMS_SIZE; i++){
        if(algorithm == algorithms[i][0]){
            valid = true;
            break;
        }
    }
    if(!valid)
        print_error("Algorithm is not available\n");
    
}

/**
 * @brief
 * Print the algorithm name
 * @param algorithm the algorithm to print
 * @return the name of the algorithm
 */
const char* print_algorithm(char algorithm){
    for(int i = 0; i < ALGORITHMS_SIZE; i++){
        if(algorithm == algorithms[i][0]){
            return algorithms[i];
        }
    }

    return "Unknown";
}

/**
 * @brief
 * Print the available algorithms
 */
void print_algorithms(){
    for(int i = 0; i < ALGORITHMS_SIZE; i++){
        printf("\t%s\n", algorithms[i]);
    }
}

/**
 * @brief
 * Print the available parameters
 */
void print_parameters(){
    for(int i = 0; i < PARAMS; i++){
        printf("\t%s\n", parameters[i]);
    }
}

void setAlgorithmId(instance *inst, char *algorithmID){
    switch(inst->algorithm){
        case 'V':
            snprintf(algorithmID, sizeof(algorithmID), "%c_%d_%d", inst->algorithm, inst->params[KICK], inst->params[K_OPT]);
            break;
        case 'G':
            snprintf(algorithmID, sizeof(algorithmID), "%c_%d_%d", inst->algorithm, inst->params[ALPHA], inst->params[MIN_COSTS]);
            break;
        case 'T':
            snprintf(algorithmID, 1000, "%c_%d_%d_%d", inst->algorithm, inst->params[MIN_TENURE], inst->params[MAX_TENURE], inst->params[TENURE_STEP]);
            break;
        case 'B':
            snprintf(algorithmID, sizeof(algorithmID), "%c_%lf", toupper(inst->algorithm), inst->time_limit); 
            break;
        default:
            print_error("Algorithm not implemented");
            break;
    }

}