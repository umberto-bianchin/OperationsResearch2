#define _CRT_SECURE_NO_WARNINGS
#include <parsers.h>

/**
 * @brief
 * Parses the command line arguments and sets the instance parameters
 */
void parse_command_line(int argc, char** argv, instance *inst) { 	
	if (VERBOSE >= DEBUG)
		printf(" running %s with %d parameters \n", argv[0], argc-1); 

    int help = 0; 
	if(argc <= 1) 
		help = 1;

	for ( int i = 1; i < argc; i++ ) 	{ 		
        if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 					// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 					// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 						// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->time_limit = atof(argv[++i]); continue; }				// total time limit
		if ( strcmp(argv[i],"-t") == 0 ) { inst->time_limit = atof(argv[++i]); continue; }						// total time limit
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->seed = abs(atoi(argv[++i])); continue; } 					// random seed
		if ( strcmp(argv[i],"-n") == 0 ) { inst->nnodes= atoi(argv[++i]); continue; } 							// max n. of nodes
		if ( strcmp(argv[i],"-nodes") == 0 ) { inst->nnodes= atoi(argv[++i]); continue; } 						// max n. of nodes
		if ( strcmp(argv[i],"-a") == 0 ) { inst->algorithm = toupper(argv[++i][0]); continue; } 				// algorithm to use
		if ( strcmp(argv[i],"-algorithm") == 0 ) { inst->algorithm = toupper(argv[++i][0]); continue; } 		// algorithm to use
		if ( strcmp(argv[i],"-r") == 0 ) { inst->running_mode = tolower(argv[++i][0]); continue; } 				// running mode
		if ( strcmp(argv[i],"-i") == 0 ) { inst->integer_costs = atoi(argv[++i]); continue; } 					// integer costs
		if ( strcmp(argv[i],"-kick") == 0 ) { inst->params[KICK] = atoi(argv[++i]); continue; } 			 	// kick param
		if ( strcmp(argv[i],"-kopt") == 0 ) { inst->params[K_OPT] = atoi(argv[++i]); continue; } 			 	// kick param
		if ( strcmp(argv[i],"-alpha") == 0 ) { inst->params[ALPHA] = atoi(argv[++i]); continue; } 				// alpha param
		if ( strcmp(argv[i],"-minc") == 0 ) { inst->params[MIN_COSTS] = atoi(argv[++i]); continue; } 			// min_costs param
		if ( strcmp(argv[i],"-maxt") == 0 ) { inst->params[MAX_TENURE] = atoi(argv[++i]); continue; } 			// max_tenure param
		if ( strcmp(argv[i],"-mint") == 0 ) { inst->params[MIN_TENURE] = atoi(argv[++i]); continue; } 			// min_tenure param
		if ( strcmp(argv[i],"-stept") == 0 ) { inst->params[TENURE_STEP] = atoi(argv[++i]); continue; } 		// tenure_step param
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 											// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 											// help
		help = 1;
    }      

	// print current parameters
	if (VERBOSE >= ERROR){	
		printf("available parameters --------------------------------------------------\n");
		printf("-file %s\n", inst->input_file); 
        printf("-time_limit %lf\n", inst->time_limit);
		printf("-seed %d\n", inst->seed); 
		printf("-n %d\n", inst->nnodes); 
		printf("-a %c\n", inst->algorithm);
		printf("-i %d\n", inst->integer_costs);
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}        
	
	if(help == 1) {
		printf("----------------------------------------------------------------------------------------------\n");
		printf("To use correctly this program you have to insert as a parameter the -file or -nodes described here\n");
		printf("Type the command with the following parameters:\n");
		printf("-file <file_name> : the path of the file containing the TSP instance\n");
		printf("-t || -time_limit <time> : the time limit in seconds for the algorithm\n");
		printf("-seed <seed> : the seed for the random number generator\n");
		printf("-n || -nodes <nodes> : the number of nodes random generated using seed\n");
		printf("-i [0, 1] : 1 if you want to use integer costs\n");
		printf("-a <algorithm> : the algorithm to use\n");
		print_algorithms();
		printf("-r <running mode> : the running mode: [B] for benchmark, [N] for normal, [C] for CPLEX\n");
		printf("-[param] <parameter> : for each algorithm you can choose different params\n");
		print_parameters();
		printf("-help : print this help\n\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
		exit(1);
	}

	check_input(inst);
}    

/**
 * @brief 
 * Checks if the input parameters are correct
 */
void check_input(instance *inst){
	if (VERBOSE >= DEBUG)
		printf("... checking input\n");

	if(inst->integer_costs != 0 && inst->integer_costs != 1)
		print_error("Integer costs can be only ON [1] or OFF [0]\n Use the -help command to see how to select an algorithm"); 

	if(inst->algorithm == ' ' && inst->running_mode != 'c')
		print_error("Algorithm not set!\n Use the -help command to see how to select an algorithm"); 
	
	if(inst->running_mode != 'c')
		check_valid_algorithm(inst->algorithm);
	
	if(inst->nnodes <= 0 && (strcmp(inst->input_file, "NULL") == 0))
		print_error("Invalid options, use the -help command to see how to run this program.");
	
	if(inst->running_mode != 'b' && inst->running_mode != 'n' && inst->running_mode != 'c')
		print_error("Invalid running mode, use the -help command to see how to run this program.");
	
	if(inst->params[K_OPT] < 3 && inst->algorithm == 'V')
		print_error("Invalid kopt choice, use the -help command to see how to run this program.");
}

/**
 * @brief
 * Reads the input file and sets the instance parameters,
 * such as the number of nodes, the coordinates of each node and the costs between each pair of nodes.
 * Is a simplified CVRP parser, not all SECTIONs are detected.
 */
void read_input(instance *inst) {    
	// if the number of nodes is already setted, the input file is not read and the coordinates are generated randomly
	if(inst->nnodes > 0){
		allocate_instance(inst);
		set_random_coord(inst);
		return;
	}

	FILE *fin = fopen(inst->input_file, "r");
	
	if(fin == NULL)
		print_error("Input file not found!");
	
	inst->nnodes = -1;

	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 
	
	int do_print = (VERBOSE >= DEBUG);

	while(fgets(line, sizeof(line), fin) != NULL){
		if (VERBOSE >= DEBUG){ printf("%s",line); fflush(NULL); }
		if (strlen(line) <= 1) continue; // skip empty lines
	    par_name = strtok(line, " :");
		if (VERBOSE >= DEBUG){ printf("parameter \"%s\" ",par_name); fflush(NULL); }

		if(strncmp(par_name, "NAME", 4) == 0){
			active_section = 0;
			continue;
		}

		if(strncmp(par_name, "COMMENT", 7) == 0){
			active_section = 0;   
			token1 = strtok(NULL, "");  
			//if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
			continue;
		}   
		
		if(strncmp(par_name, "TYPE", 4) == 0){
			token1 = strtok(NULL, " :");  
			if ( strncmp(token1, "TSP",3) != 0 ) print_error(" format error:  only TYPE == CVRP implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "DIMENSION", 9) == 0 ) {
			if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes); 
			//inst->demand = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			allocate_instance(inst);
			active_section = 0;  
			continue;
		}         
		
		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ){
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;   
			continue;
		}

		if ( strncmp(par_name, "EOF", 3) == 0 ){
			active_section = 0;
			break;
		}
		
			
		if ( active_section == 1 ){ // within NODE_COORD_SECTION
			int i = atoi(par_name) - 1; 
			if ( i < 0 || i >= inst->nnodes ) print_error(" ... unknown node in NODE_COORD_SECTION section");     
			token1 = strtok(NULL, " :,");
			token2 = strtok(NULL, " :,");
			inst->xcoord[i] = atof(token1);
			inst->ycoord[i] = atof(token2);
			if ( do_print ) printf(" ... node %4d at coordinates ( %15.7lf , %15.7lf )\n", i+1, inst->xcoord[i], inst->ycoord[i]); 
			continue;
		}       
	}                

	fclose(fin);    
	
}