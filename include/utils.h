#ifndef UTILS_H_  
#define UTILS_H_

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#define INFO		20
#define ERROR		50
#define DEBUG		100
#define VERBOSE		INFO
#define MAX_NO_IMPROVEMENT 10000

void print_error(const char *err, bool terminate);

#endif /* UTILS_H_ */ 