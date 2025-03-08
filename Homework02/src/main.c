#include <tsp_utilities.h>
#include <heuristics.h>
#include <parsers.h>
#include <chrono.h>

int main(int argc, char **argv) {
	char algorithm;
	double t1, t2;
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	instance inst;

	parse_command_line(argc, argv, &inst);
	  
	read_input(&inst);

	printf("Choose the algorithm to use: N for nearest neighbour, E for extra-mileage\n");
	algorithm = getchar();
	getchar();

	compute_all_costs(&inst);

	if(algorithm == 'N')
	{	
		printf("Solving problem with nearest neighbour algorithm\n");

		t1 = second(); 
		all_nearest_neighbours(&inst);
		t2 = second();
	} else if (algorithm == 'E'){
		printf("Solving problem with extra mileage algorithm\n");

		t1 = second();
		extra_mileage(&inst);
		t2 = second();
	} else{
		printf("Algorithm %c is not available\n", algorithm);
		exit(EXIT_FAILURE);
	}


	plot_solution(&inst, 1);
	
	if ( VERBOSE >= 1 ){ printf("... VRP problem solved in %lf sec.s\n", t2-t1); }
	
	free_instance(&inst);

	return 0;
}