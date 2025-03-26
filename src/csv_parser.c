#define _CRT_SECURE_NO_WARNINGS
#include <csv_parser.h>

/**
 * @brief 
 * Parse a csv file and store the data in a 3D array
 * @param csvData the 3D array
 */
int parse_csv(char ***csvData)
{
    FILE *fp = fopen(FILENAME, "r");
    if (!fp) {
        return 0;
    }

    char line[MAX_LINE_LEN];
    int row = 0, numCols = 0;
    
    while (fgets(line, MAX_LINE_LEN, fp) != NULL && row < MAX_ROWS) {
        line[strcspn(line, "\n")] = '\0';

        int col = 0;
        char *token = strtok(line, ",\n");
        token = strtok(NULL, ",\n");        //skip first column
        while (token && col < MAX_COLS) {
            strncpy(csvData[row][col], token, MAX_LINE_LEN);
            csvData[row][col][MAX_LINE_LEN - 1] = '\0';
            col++;
            token = strtok(NULL, ",\n");
        }

        if (row == 0) {
            numCols = col;
        }

        // Clear the rest of columns in this row
        for (; col < MAX_COLS; col++) {
            csvData[row][col][0] = '\0';
        }

        row++;
    }

    fclose(fp);
    
    return numCols;
}

/**
 * @brief 
 * Call the function to read the csv file, then write in it the new data
 * @param bestCosts the new costs to insert in the csv file
 * @param algorithmID the algorithm ID for the header
 */
void write_csv(double bestCosts[MAX_ROWS - 1], char *algorithmID)
{   
    // 3D array for CSV: csvData[row][col][string]
    char ***csvData = malloc(MAX_ROWS * sizeof(char **));
    for (int r = 0; r < MAX_ROWS; r++) {
        csvData[r] = malloc(MAX_COLS * sizeof(char *));
        for (int c = 0; c < MAX_COLS; c++) {
            csvData[r][c] = calloc(MAX_LINE_LEN, sizeof(char));
            csvData[r][c][0] = '\0';
        }
    }

    // read existing data to append new data in the same file
    int cols = parse_csv(csvData);

    strncpy(csvData[0][cols], algorithmID, MAX_LINE_LEN);

    for (int i = 0; i < MAX_ROWS - 1; i++) {
        snprintf(csvData[i + 1][cols], MAX_LINE_LEN, "%.6f", bestCosts[i]);
    }

    FILE *fp = fopen(FILENAME, "w");
    if (!fp) {
        printf("Cannot open CSV for writing\n");
        return;
    }

    for (int r = 0; r < MAX_ROWS; r++) {
        if (csvData[r][0][0] == '\0') {
            break;
        }
        
        if(r == 0)
            fprintf(fp, "%d,", cols+1);
        else
            fprintf(fp, "tsp_instance%d,", r);

        // For each column in this row until we see an empty cell
        for (int c = 0; c < cols + 1; c++) {
            if (csvData[r][c][0] == '\0') break;
            fprintf(fp, "%s,", csvData[r][c]);
        }
        fprintf(fp, "\n");
    }

    fclose(fp);
    for (int r = 0; r < MAX_ROWS; r++) {
        for (int c = 0; c < MAX_COLS; c++) {
            free(csvData[r][c]);
        }
        free(csvData[r]);
    }
    free(csvData);
    
    printf("Results are in file '%s', column %d\n", FILENAME, cols + 1);
}