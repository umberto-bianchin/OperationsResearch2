#include "tsp_utilities.h"

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }

void free_instance(instance *inst)
{     
	free(inst->xcoord);
	free(inst->ycoord);
}

void choose_rand_sol(instance *inst)
{
    inst->best_sol = (double *)malloc(5 * sizeof(double));

    srand(inst->seed);
    for (int i = 0; i < 5; i++) {
        int index = rand() % inst->nnodes;
        inst->best_sol[i] = index;
    }

    if(VERBOSE >= 20) 
    {
        printf("Choosen value for best_sol: ");
        for (int i = 0; i < 5; i++) printf("%f ", inst->best_sol[i]);
		printf("\n");
    }
}

void plot_solution(instance *inst)
{
    #ifdef WIN32
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
    fprintf(gnuplotPipe, "plot '-' with linespoints linestyle 1 linewidth 2 pointtype 7 pointsize 1.5 linecolor 'blue'\n");

	for(int i=0; i<5; i++)
	{
		int idx = inst->best_sol[i];
		fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[idx], inst->ycoord[idx]);
	}
    
	fprintf(gnuplotPipe, "%lf %lf\n", inst->xcoord[(int)inst->best_sol[0]], inst->ycoord[(int)inst->best_sol[0]]);

    fprintf(gnuplotPipe, "e\n");

    #ifdef WIN32
		_pclose(gnuplotPipe);
	#else
		pclose(gnuplotPipe);
	#endif

}