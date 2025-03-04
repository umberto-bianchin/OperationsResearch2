#include "../include/tsp_utilities.h"

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }

void free_instance(instance *inst)
{     
	free(inst->xcoord);
	free(inst->ycoord);
}

void plot_solution(instance *inst, int *solution)
{
    #ifdef _WIN32
		FILE *gnuplotPipe = _popen("gnuplot -persistent", "w");
	#else
		FILE *gnuplotPipe = popen("gnuplot -persistent", "w");
	#endif

	fprintf(gnuplotPipe, "set title 'TSP Solution'\n");
    fprintf(gnuplotPipe, "set xlabel 'X'\n");
    fprintf(gnuplotPipe, "set ylabel 'Y'\n");
    fprintf(gnuplotPipe, "set grid\n");
	fprintf(gnuplotPipe, "set key outside top\n");
	fprintf(gnuplotPipe, "set term qt title 'TSP Solution'\n");
    fprintf(gnuplotPipe, "plot '-' with linespoints linestyle 1 linewidth 2 pointtype 7 pointsize 1.5 linecolor 'blue' title 'TSP Solution'\n");

	for(int i = 0; i < inst->nnodes; i++)
	{
		int idx = solution[i];
		fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
	}
    
	fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[(int)solution[0]], inst->ycoord[(int)solution[0]]);

    fprintf(gnuplotPipe, "e\n");

    #ifdef _WIN32
		_pclose(gnuplotPipe);
	#else
		pclose(gnuplotPipe);
	#endif
}

bool check_solution(int sol, double *cost, instance *inst)
{
	return false;
}

void update_best_sol(int sol, double *cost, instance *inst)
{
	
}

double dist(int i, int j, instance *inst)
{
	return 0.0;
}