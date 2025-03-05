#include <tsp_utilities.h>
#include <parsers.h>
#include <chrono.h>

int main(int argc, char **argv) 
{ 
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	double t1 = second(); 
	instance inst;

	parse_command_line(argc, argv, &inst);     
	
	//printf(" file %s has %d non-empty lines\n", inst.input_file, number_of_nonempty_lines(inst.input_file)); exit(1);
	  
	read_input(&inst);  
	//if ( VRPopt(&inst) ) print_error(" error within VRPopt()");
    int *random_solution = (int *)malloc(5 * sizeof(int));
	choose_rand_sol(&inst);
	double t2 = second();

	check_solution(&inst);
	plot_solution(&inst);
	
	if ( VERBOSE >= 1 ){ printf("... VRP problem solved in %lf sec.s\n", t2-t1); }
	
	free_instance(&inst);

	return 0;
}