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
	  
	read_input(&inst);

	compute_all_costs(&inst);

	//choose_rand_sol(&inst);
	for(int i = 0; i < inst.nnodes; i++)
		inst.solution[i] = i;

		inst.solution[inst.nnodes] = 0;

	calc_solution_cost(&inst);
	two_opt(&inst);
	double t2 = second();

	check_solution(&inst, 0);
	plot_solution(&inst, 0);
	
	if ( VERBOSE >= 1 ){ printf("... VRP problem solved in %lf sec.s\n", t2-t1); }
	
	free_instance(&inst);

	return 0;
}