#define _CRT_SECURE_NO_WARNINGS
#include <utils.h>
#include <parsers.h>
#include <cplex_utilities.h>
#include <matheuristics.h>

/**
 * @brief Print an error message and terminate the program.
 *
 * @param err Null-terminated error string to display.
 */
void print_error(const char *err){
    printf("\n\n ERROR: %s \n\n", err); 
    fflush(NULL);
    
    exit(EXIT_FAILURE);
    
}


/**
 * @brief Initialize a dynamic solutions history struct.
 *
 * Allocates arrays for costs and iteration times with an initial capacity.
 *
 * @param sol Pointer to solutions struct to set up.
 */
void allocate_solution_struct(solutions *sol){
    sol->capacity = 16;
    sol->size = 0;
    sol->all_costs = (double*)malloc(sol->capacity * sizeof(double));
    sol->iteration_times = (double*)malloc(sol->capacity * sizeof(double));
}

/**
 * @brief Release memory held by a solutions history struct.
 *
 * Frees both the cost and time arrays if non-NULL.
 *
 * @param sol Pointer to solutions struct to clean up.
 */
void free_solution_struct(solutions *sol){
    if(sol->all_costs != NULL)
        free(sol->all_costs);
    if(sol->iteration_times != NULL)
        free(sol->iteration_times);
}

/**
 * @brief Append a new solution record, expanding storage if needed.
 *
 * If size reaches capacity, arrays grow by 20 elements via realloc.
 *
 * @param sol   Solutions struct to update.
 * @param cost  Cost value to append.
 * @param time  Timestamp (or -1 to skip recording time).
 */
void add_solution(solutions *sol, double cost, double time){
    if(sol->size == sol->capacity){
        sol->capacity += 20;
        // expand arrays, preserving existing data
        sol->all_costs = (double*)realloc(sol->all_costs, sol->capacity * sizeof(double));
        sol->iteration_times = (double*)realloc(sol->iteration_times, sol->capacity * sizeof(double));
    }

    sol->all_costs[sol->size] = cost;
    if(time != -1)
        sol->iteration_times[sol->size] = time;
    
    sol->size++;
}

/**
 * @brief Plot a single TSP tour using gnuplot, saving to PNG file.
 *
 * Builds filename with algorithm ID and timestamp, then streams commands
 * to gnuplot to draw edges and nodes.
 *
 * @param inst Pointer to TSP instance with coords and best_solution.
 * @param s    Solution path to render (must have path array).
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

/**
 * @brief Plot cost histories: best vs. all over iterations.
 *
 * Saves a PNG showing two line series for solution tracking.
 *
 * @param inst Pointer to TSP instance with history_best_costs and history_costs.
 */
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

/**
 * @brief Plot CPLEX solution costs over time.
 *
 * Uses history_best_costs.iteration_times for X-axis.
 *
 * @param inst Pointer to TSP instance with CPLEX history data.
 */
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
     fprintf(gnuplotPipe, "set title 'Algorithm: %s, Solution Cost: %.4lf, Time consumed: %.2lf'\n", print_algorithm(inst->algorithm), inst->best_solution.cost, (second() - inst->t_start));
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
 * @brief Select and execute the TSP algorithm indicated in the instance.
 *
 * Prompts for time limit if not set for certain methods, then plots results.
 *
 * @param inst Pointer to initialized TSP instance (algorithm, time_limit set).
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
        two_opt(inst, &s, inst->time_limit);        // warm-up 2-opt

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
        case 'C':
            if(inst->time_limit == INF_COST){
                printf("Please, insert time limit (seconds): ");
                scanf("%lf", &inst->time_limit);
                getchar();
                printf("Updated time to solve this problem: %lf seconds\n", inst->time_limit);
            }
            if(inst->algorithm == 'B')
                printf("Solving problem with benders\n");
            else
                printf("Solving problem with branch and cut\n");
            TSPopt(inst);
            break;
        case 'H':
            printf("Solving problem with hard fixing\n");
            hard_fixing(inst);
            break;
        default:
            print_error("Algorithm is not available\n");
            break;
    }

	t2 = second();
	if (VERBOSE >= INFO)
		printf("\nTSP problem solved in %lf sec.s\n", t2-inst->t_start);
    
	plot_solution(inst, &(inst->best_solution));

    if(inst->algorithm != 'B' && inst->algorithm != 'C'){
        plot_solutions(inst);
    } else{
        plot_cplex_solutions(inst);
    }
}

/**
 * @brief Benchmark CPLEX-based algorithm over random instances, recording times.
 *
 * Generates MAX_ROWS-1 random TSPs, solves with TSPopt, and writes CSV.
 *
 * @param inst Pointer to instance (algorithm must be 'B' or 'C').
 */
void benchmark_algorithm_by_time(instance *inst){
    double solvingTimes[MAX_ROWS - 1];

    for (int i = 0; i < MAX_ROWS - 1; i++) {
        inst->seed = i * 1000;
        inst->best_solution.cost = INF_COST;

        set_random_coord(inst);
        compute_all_costs(inst);
        inst->t_start = second();

        if(inst->algorithm != 'B' && inst->algorithm != 'C'){
            print_error("Impossible benchmarking algorithm different from Benders or Branch and Cut by time");
        }

        printf("\nRunning CPLEX TSP solver on random coordinates with seed %d...\n", inst->seed);        

        TSPopt(inst);

        double elapsed_time = second() - inst->t_start;
        solvingTimes[i] = elapsed_time;
    }

    char algorithmID[64];
    setAlgorithmId(inst, algorithmID);

    write_csv(solvingTimes, algorithmID, inst->algorithm); 
}

/**
 * @brief Benchmark heuristic algorithms, recording best costs over seeds.
 *
 * Runs specified method MAX_ROWS-1 times with varied seeds, then CSV.
 *
 * @param inst Pointer to instance with algorithm and time_limit.
 */
void benchmark_algorithm_by_params(instance *inst){   
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
            case 'H':
                printf("Running hard fixing on random coordinates with seed %d...\n", inst->seed);
                hard_fixing(inst);
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
 * @brief Validate a chosen algorithm character against available list.
 *
 * @param algorithm Character code to check (e.g., 'N','E', etc.).
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
 * @brief Return human-readable name for an algorithm code.
 *
 * @param algorithm Single-letter code.
 * @return Pointer to static string name or "Unknown".
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
 * @brief Print all available algorithm names to stdout.
 */
void print_algorithms(){
    for(int i = 0; i < ALGORITHMS_SIZE; i++){
        printf("\t%s\n", algorithms[i]);
    }
}

/**
 * @brief Print all configurable parameter names to stdout.
 */
void print_parameters(){
    for(int i = 0; i < PARAMS; i++){
        printf("\t%s\n", parameters[i]);
    }
}

/**
 * @brief Build a compact run identifier string based on algorithm and params.
 *
 * Fills algorithmID buffer with code_param1_param2... pattern.
 *
 * @param inst        TSP instance with algorithm char and params set.
 * @param algorithmID Output buffer (char[]) to write identifier into.
 */
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
            snprintf(algorithmID, 1000, "%c_%.2lf_%d_%d", toupper(inst->algorithm), inst->time_limit, inst->params[WARMUP], inst->params[POSTING]); 
            break;
        case 'C':
            snprintf(algorithmID, 1000, "%c_%.2lf_%d_%d_%d_%d", toupper(inst->algorithm), inst->time_limit, inst->params[WARMUP], inst->params[POSTING], inst->params[DEPTH], inst->params[CONCORDE]); 
            break;
        case 'H':
            snprintf(algorithmID, 1000, "%c_%.2lf_%d_%d", toupper(inst->algorithm), inst->time_limit, inst->params[FIXEDPROB], inst->params[PROBABILITY]); 
            break;           
        default:
            print_error("Algorithm not implemented");
            break;
    }

}