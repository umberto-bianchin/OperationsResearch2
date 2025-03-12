#include <tsp_utilities.h>
#include <parsers.h>
#include <chrono.h>

void debug5nodes(){
	instance inst;

	initialize_instance(&inst);
	compute_all_costs(&inst);

	all_nearest_neighbours(inst);

	for(int i = 0; i < inst.nnodes + 1; i++)
		printf("%d ", inst.solution[i]);
	printf("\n");

	plot_solution(&inst, 1);
	
	free_instance(&inst);
}

int main(int argc, char **argv) {
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(EXIT_FAILURE); }       

	instance inst;

	initialize_instance(&inst);
	
	/* Debug code */
	debug5nodes();
	
	/* normal code 
	parse_command_line(argc, argv, &inst);
	read_input(&inst);

	compute_all_costs(&inst);

	choose_run_algorithm(&inst);

	plot_solution(&inst, 1);
	
	free_instance(&inst);
	*/
	return 0;
}