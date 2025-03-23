#include <parsers.h>
#include <chrono.h>

#include <heuristics.h>
#include <tsp_utilities.h>
#include <utils.h>

int main(int argc, char **argv) {

    printf("Choose run mode:\n");
    printf(" - 'b' for benchmark by params\n");
    printf(" - 'n' for normal execution\n");

    char mode = tolower(getchar());
    getchar();

    if (mode == 'b') {
        printf("Running in BENCHMARK mode.\n\n");
        benchmark_algorithm_by_params();
    } else if (mode == 'n') {
        printf("Running in NORMAL mode.\n\n");
        instance inst;

	    initialize_instance(&inst);

	    parse_command_line(argc, argv, &inst);
	    read_input(&inst);

	    compute_all_costs(&inst);
        choose_run_algorithm(&inst);
        free_instance(&inst);
    } else {
        printf("Unknown mode '%c'. Exiting.\n", mode);
        exit(EXIT_FAILURE);
    }
	
	return 0;
}