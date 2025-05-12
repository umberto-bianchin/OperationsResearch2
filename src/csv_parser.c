#define _CRT_SECURE_NO_WARNINGS
#include <csv_parser.h>
#include <utils.h>

/**
 * @brief
 * Parse a CSV file (skipping its first column) into a preallocated 3D char array.
 *
 * Reads up to MAX_ROWS lines and MAX_COLS columns per line. Each token
 * is copied into csvData[row][col] and null-terminated.
 *
 * @param csvData  Preallocated 3D array: csvData[row][col][char]
 * @param fileName Path to the CSV file to read.
 * @return Number of columns parsed (excluding the skipped first column),
 *         or 0 if the file could not be opened.
 */
int parse_csv(char ***csvData, char *fileName){ 
    FILE *fp = fopen(fileName, "r");
    if (!fp) {
        return 0;
    }

    char line[MAX_LINE_LEN];
    int row = 0, numCols = 0;
    
    while (fgets(line, MAX_LINE_LEN, fp) != NULL && row < MAX_ROWS) {
        line[strcspn(line, "\n")] = '\0';

        int col = 0;
        // tokenize by comma; skip first cell
        char *token = strtok(line, ",\n");
        token = strtok(NULL, ",\n");

        while (token && col < MAX_COLS) {
            strncpy(csvData[row][col], token, MAX_LINE_LEN);
            csvData[row][col][MAX_LINE_LEN - 1] = '\0';
            col++;
            token = strtok(NULL, ",\n");
        }

        if (row == 0) {
            numCols = col;
        }

        // clear any unused columns in this row
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
 * Append a new column of data (bestCosts) to an existing CSV file.
 *
 * Calls parse_csv() to load existing contents into memory, then writes
 * back with algorithmID in header and formatted costs in new column.
 *
 * @param bestCosts     Array of double values, length MAX_ROWS-1, to append.
 * @param algorithmID   String ID to place in header for new column.
 * @param alg           Single-character algorithm code (for file naming).
 */
void write_csv(double bestCosts[MAX_ROWS - 1], char *algorithmID, char alg){   
    char fileName[128];

    // choose file name template based on configuration
    if(MULTIPLE_FILES){
        snprintf(fileName, sizeof(fileName), "%s_%c.csv", FILENAME, alg);
    } else{
        snprintf(fileName, sizeof(fileName), "%s.csv", FILENAME);
    }

    /// allocate 3D array: rows x cols x string
    char ***csvData = malloc(MAX_ROWS * sizeof(char **));
    for (int r = 0; r < MAX_ROWS; r++) {
        csvData[r] = malloc(MAX_COLS * sizeof(char *));
        for (int c = 0; c < MAX_COLS; c++) {
            // each cell can hold up to MAX_LINE_LEN chars
            csvData[r][c] = calloc(MAX_LINE_LEN, sizeof(char));
            csvData[r][c][0] = '\0';
        }
    }

    // read existing CSV into csvData
    int cols = parse_csv(csvData, fileName);

    // place algorithmID in header at new column index
    strncpy(csvData[0][cols], algorithmID, MAX_LINE_LEN);

    // fill in bestCosts for each subsequent row
    for (int i = 0; i < MAX_ROWS - 1; i++) {
        snprintf(csvData[i + 1][cols], MAX_LINE_LEN, "%.6f", bestCosts[i]);
    }

    // open file for overwrite
    FILE *fp = fopen(fileName, "w");
    if (!fp) {
        print_error("Cannot open CSV for writing\n");
    }

    // write out each populated row
    for (int r = 0; r < MAX_ROWS; r++) {
        // if first column empty, assume end of data
        if (csvData[r][0][0] == '\0') {
            break;
        }
        
        // write row label
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
    
    printf("Results are in file '%s', column %d\n", fileName, cols + 1);
}