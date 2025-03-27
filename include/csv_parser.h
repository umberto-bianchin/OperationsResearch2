#ifndef CSV_PARSER_H_  
#define CSV_PARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ROWS 11          //NUM_FILES + 1 (header)
#define MAX_COLS 50             //maximum number of columns in the file
#define MAX_LINE_LEN 1024
#define FILENAME "utils/results"

int parse_csv(char ***csvData, char *fileName);
void write_csv(double bestCosts[MAX_ROWS - 1], char *algorithmID, char alg);

#endif   /* CSV_PARSER_H_ */ 