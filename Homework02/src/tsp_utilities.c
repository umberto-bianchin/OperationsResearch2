#include "tsp_utilities.h"

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }

void free_instance(instance *inst)
{     
	free(inst->xcoord);
	free(inst->ycoord);
	free(inst->costs);
	free(inst->solution);
	free(inst->best_sol);
}

void choose_rand_sol(instance *inst)
{
    //inst->solution = (int *)malloc(5 * sizeof(int));

    srand(inst->seed);
    for (int i = 0; i < 5; i++) {
        inst->solution[i] = rand() % inst->nnodes;
    }

    if(VERBOSE >= 20) 
    {
        printf("Choosen value for best_sol: ");
        for (int i = 0; i < 5; i++) printf("%d ", inst->solution[i]);
		printf("\n");
    }
}

void plot_solution(instance *inst, bool best)
{
    #ifdef _WIN32
		FILE *gnuplotPipe = _popen("gnuplot -persistent", "w");
	#else
		FILE *gnuplotPipe = popen("gnuplot", "w");
	#endif

	int *solution = best ? inst->best_sol : inst->solution;

	if(solution == NULL){printf("Solution is not initialized"); exit(1);}

	fprintf(gnuplotPipe, "set title 'TSP Solution'\n");
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set grid\n");
	fprintf(gnuplotPipe, "set key outside top\n");
	fprintf(gnuplotPipe, "set term qt title 'TSP Solution'\n");
    fprintf(gnuplotPipe, "plot '-' with linespoints linestyle 1 linewidth 2 pointtype 7 pointsize 1.5 linecolor 'blue' title 'TSP Solution'\n");

	for(int i = 0; i < inst->nnodes; i++)
	{
		int idx = inst->solution[i];
		fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
	}
    
	fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[(int)inst->solution[0]], inst->ycoord[(int)inst->solution[0]]);

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

bool check_solution(instance *inst)
{
	return false;
}

void update_best_sol(instance *inst)
{
	
}

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