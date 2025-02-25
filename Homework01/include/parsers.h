#ifndef PARSERS_H_  

#define PARSERS_H_

#include <tsp.h>

void parse_command_line(int argc, char** argv, instance *inst);
void read_input(instance *inst);
void print_error(const char *err);

#endif   /*PARSERS_H_ */ 