#include <parsers.h>

void parse_command_line(int argc, char** argv, instance *inst) 
{ 	
	if ( VERBOSE >= 100 ) printf(" running %s with %d parameters \n", argv[0], argc-1); 
		
	// default   
	inst->seed = 0;
	inst->timelimit = CPX_INFBOUND;
	inst->nnodes = -1;
	strcpy(inst->input_file, "NULL");

    int help = 0; if ( argc < 1 ) help = 1;	
	for ( int i = 1; i < argc; i++ ) 
	{ 		
        if ( strcmp(argv[i],"-file") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-input") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 			// input file
		if ( strcmp(argv[i],"-f") == 0 ) { strcpy(inst->input_file,argv[++i]); continue; } 				// input file
		if ( strcmp(argv[i],"-time_limit") == 0 ) { inst->timelimit = atof(argv[++i]); continue; }		// total time limit
		if ( strcmp(argv[i],"-seed") == 0 ) { inst->seed = abs(atoi(argv[++i])); continue; } 		// random seed
		if ( strcmp(argv[i],"-n") == 0 ) { inst->nnodes= atoi(argv[++i]); continue; } 					// max n. of nodes
		if ( strcmp(argv[i],"-nodes") == 0 ) { inst->nnodes= atoi(argv[++i]); continue; } 					// max n. of nodes
		if ( strcmp(argv[i],"-help") == 0 ) { help = 1; continue; } 									// help
		if ( strcmp(argv[i],"--help") == 0 ) { help = 1; continue; } 									// help
		help = 1;
    }      

	if ( help || (VERBOSE >= 10) )		// print current parameters
	{
		printf("\n\navailable parameters (vers. 16-may-2015) --------------------------------------------------\n");
		printf("-file %s\n", inst->input_file); 
        printf("-time_limit %lf\n", inst->timelimit);
		printf("-seed %d\n", inst->seed); 
		printf("-n %d\n", inst->nnodes); 
		printf("\nenter -help or --help for help\n");
		printf("----------------------------------------------------------------------------------------------\n\n");
	}        
	
	if ( help ) exit(1);

}    

void read_input(instance *inst) // simplified CVRP parser, not all SECTIONs detected  
{
                            
	FILE *fin = fopen(inst->input_file, "r");
	if ( fin == NULL ) print_error(" input file not found!");
	
	inst->nnodes = -1;

	char line[180];
	char *par_name;   
	char *token1;
	char *token2;
	
	int active_section = 0; // =1 NODE_COORD_SECTION, =2 DEMAND_SECTION, =3 DEPOT_SECTION 
	
	int do_print = ( VERBOSE >= 1000 );

	while ( fgets(line, sizeof(line), fin) != NULL ) 
	{
		if ( VERBOSE >= 2000 ) { printf("%s",line); fflush(NULL); }
		if ( strlen(line) <= 1 ) continue; // skip empty lines
	    par_name = strtok(line, " :");
		if ( VERBOSE >= 3000 ) { printf("parameter \"%s\" ",par_name); fflush(NULL); }

		if ( strncmp(par_name, "NAME", 4) == 0 ) 
		{
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "COMMENT", 7) == 0 ) 
		{
			active_section = 0;   
			token1 = strtok(NULL, "");  
			//if ( VERBOSE >= 10 ) printf(" ... solving instance %s with model %d\n\n", token1, inst->model_type);
			continue;
		}   
		
		if ( strncmp(par_name, "TYPE", 4) == 0 ) 
		{
			token1 = strtok(NULL, " :");  
			if ( strncmp(token1, "TSP",3) != 0 ) print_error(" format error:  only TYPE == CVRP implemented so far!!!!!!"); 
			active_section = 0;
			continue;
		}

		if ( strncmp(par_name, "DIMENSION", 9) == 0 ) 
		{
			if ( inst->nnodes >= 0 ) print_error(" repeated DIMENSION section in input file");
			token1 = strtok(NULL, " :");
			inst->nnodes = atoi(token1);
			if ( do_print ) printf(" ... nnodes %d\n", inst->nnodes); 
			//inst->demand = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->xcoord = (double *) calloc(inst->nnodes, sizeof(double)); 	 
			inst->ycoord = (double *) calloc(inst->nnodes, sizeof(double));
			inst->best_sol = (int *) calloc(inst->nnodes + 1, sizeof(int));
			inst->solution = (int *) calloc(inst->nnodes + 1, sizeof(int));  
			inst->costs = (double *) calloc(inst->nnodes, sizeof(double));
			
			active_section = 0;  
			continue;
		}         
		
		if ( strncmp(par_name, "NODE_COORD_SECTION", 18) == 0 ) 
		{
			if ( inst->nnodes <= 0 ) print_error(" ... DIMENSION section should appear before NODE_COORD_SECTION section");
			active_section = 1;   
			continue;
		}

		if ( strncmp(par_name, "EOF", 3) == 0 ) 
		{
			active_section = 0;
			break;
		}
		
			
		if ( active_section == 1 ) // within NODE_COORD_SECTION
		{
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