#include <utils.h>

/**
 * @brief
 * Prints an error message
 * @param terminate if = 1 the program terminates its execution
 */
void print_error(const char *err, bool terminate){
    printf("\n\n ERROR: %s \n\n", err); 
    fflush(NULL);
    
    if(terminate){
        exit(EXIT_FAILURE);
    }
}