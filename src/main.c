#include <parsers.h>
#include <chrono.h>

#include <heuristics.h>
#include <tsp_utilities.h>
#include <utils.h>

int main(int argc, char **argv) {
    instance inst;

    initialize_instance(&inst);

    parse_command_line(argc, argv, &inst);

    if (inst.running_mode == 'b') {
        printf("Running in BENCHMARK mode.\n\n");
        allocate_instance(&inst);
        benchmark_algorithm_by_params(&inst);
    } else if (inst.running_mode == 'n') {
        printf("Running in NORMAL mode.\n\n");       
	    read_input(&inst);
	    compute_all_costs(&inst);
        choose_run_algorithm(&inst);
    } else {
        printf("Unknown mode '%c'. Exiting.\n", inst.running_mode);
        exit(EXIT_FAILURE);
    }
    
    free_instance(&inst);
	
	return 0;
}