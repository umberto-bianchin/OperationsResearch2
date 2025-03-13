#include <tsp_utilities.h>
#include <parsers.h>
#include <chrono.h>

#include <heuristics.h>
#include <tsp_utilities.h>

int main(int argc, char **argv) {
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(EXIT_FAILURE); }       

	instance inst;

	initialize_instance(&inst);
	
	parse_command_line(argc, argv, &inst);
	read_input(&inst);

	compute_all_costs(&inst);

	choose_run_algorithm(&inst);

	plot_solution(&inst, 1);
	
	free_instance(&inst);
	
	return 0;
}