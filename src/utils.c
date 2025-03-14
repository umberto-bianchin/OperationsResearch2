#include <utils.h>

/**
 * @brief
 * Prints an error message
 * @param err the error string
 * @param terminate true if we want the program to terminate its execution
 */
void print_error(const char *err, bool terminate){
    printf("\n\n ERROR: %s \n\n", err); 
    fflush(NULL);
    
    if(terminate){
        exit(EXIT_FAILURE);
    }
}