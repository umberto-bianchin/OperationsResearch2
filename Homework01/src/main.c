#include <tsp.h>
#include <parsers.h>

void read_input(instance *inst);
void parse_command_line(int argc, char** argv, instance *inst);
void generate_rand_sol(instance *inst);

void free_instance(instance *inst)
{     
	free(inst->xcoord);
	free(inst->ycoord);
}

int main(int argc, char **argv) 
{ 
	if ( argc < 2 ) { printf("Usage: %s -help for help\n", argv[0]); exit(1); }       
	if ( VERBOSE >= 2 ) { for (int a = 0; a < argc; a++) printf("%s ", argv[a]); printf("\n"); }

	double t1 = second(); 
	instance inst;

	parse_command_line(argc,argv, &inst);     
	
	//printf(" file %s has %d non-empty lines\n", inst.input_file, number_of_nonempty_lines(inst.input_file)); exit(1);
	  
	read_input(&inst);  
	//if ( VRPopt(&inst) ) print_error(" error within VRPopt()");
    generate_rand_sol(&inst);
	double t2 = second(); 
    
	if ( VERBOSE >= 1 )   
	{
		printf("... VRP problem solved in %lf sec.s\n", t2-t1);  
	}
	
	free_instance(&inst);
	return 0; 
}

void generate_rand_sol(instance *inst)
{
    inst->best_sol = (double *)malloc(5 * sizeof(double));

    srand(inst->randomseed); // Inizializza il generatore di numeri casuali
    for (int i = 0; i < 5; i++) {
        int index = rand() % inst->nnodes; // Sceglie un indice casuale
        inst->best_sol[i] = index;
    }

    if(VERBOSE >= 20) 
    {
        printf("Valori scelti per best_sol: ");
        for (int i = 0; i < 5; i++) printf("%f ", inst->best_sol[i]);
    }
}