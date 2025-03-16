#include <parsers.h>
#include <chrono.h>

#include <heuristics.h>
#include <tsp_utilities.h>
#include <utils.h>

int main(int argc, char **argv) {
	if ( argc < 2 ) {
		printf("Usage: %s -help for help\n", argv[0]);
		exit(EXIT_FAILURE);
	}       

	instance inst;

	initialize_instance(&inst);
	
	parse_command_line(argc, argv, &inst);
	read_input(&inst);

	compute_all_costs(&inst);

	printf("Choose run mode:\n");
    printf(" - 'b' for benchmark\n");
    printf(" - 'n' for normal execution\n");

    char mode = tolower(getchar());
    getchar();

	if (mode == 'b') {
        printf("Running in BENCHMARK mode.\n\n");
        benchmark_algorithm_by_time(&inst);
    } else if (mode == 'n') {
        printf("Running in NORMAL mode.\n\n");
        choose_run_algorithm(&inst);
    } else {
        printf("Unknown mode '%c'. Exiting.\n", mode);
        free_instance(&inst);
        exit(EXIT_FAILURE);
    }
	
	free_instance(&inst);
	
	return 0;
}