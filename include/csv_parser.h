#ifndef CSV_PARSER_H_  
#define CSV_PARSER_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_ROWS 11          //NUM_FILES + 1 (header)
#define MAX_COLS 50             //maximum number of columns in the file
#define MAX_LINE_LEN 1024
#define NUM_FILES 2             //number of files to run the algorithm on
#define FILENAME "utils/results.csv"

static const char *staticFileNames[NUM_FILES] = {
    "files/att48.tsp.txt",
    "files/att532.tsp.txt"
};

int parse_csv(char ***csvData);
void write_csv(double bestCosts[MAX_ROWS - 1], char *algorithmID);

#endif   /* CSV_PARSER_H_ */ 