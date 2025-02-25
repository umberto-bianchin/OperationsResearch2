#include "tsp.h"
#include <time.h>

void print_error(const char *err);
double second();

void print_error(const char *err) { printf("\n\n ERROR: %s \n\n", err); fflush(NULL); exit(1); }