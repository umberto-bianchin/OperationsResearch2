#include <tsp_utilities.h>
#include <parsers.h>
#include <chrono.h>

int main(int argc, char **argv) {
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	instance inst;

	parse_command_line(argc, argv, &inst);
	  
	read_input(&inst);

	compute_all_costs(&inst);

	choose_run_algorithm(&inst);

	plot_solution(&inst, 1);
	
	free_instance(&inst);

	return 0;
}