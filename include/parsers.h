#ifndef PARSERS_H_  
#define PARSERS_H_

#include <tsp_utilities.h>
#include <utils.h>

void parse_command_line(int argc, char** argv, instance *inst);
void read_input(instance *inst);
void check_input(instance *inst);

#endif   /*PARSERS_H_ */ 