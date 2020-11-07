#include <stdio.h>

/**
 * Parse a matrix
 *
 * @param fp The file descriptor to load the information from
 * @param matrixSize The size of the matrix
 * @param matrixResult Where to put the parsed matrix results (This should already be initialized to the needed size)
 * @return Whether the method was successful
 */
int parseMatrix(FILE *fp, unsigned int matrixSize, unsigned int (*matrix)[matrixSize]);

/**
 * Print a matrix
 *
 * @param fp The file descriptor to save the information to
 * @param matrixSize The size of the matrix
 * @param matrix The matrix to print
 * @return Whether the method wa successful
 */
int printMatrix(FILE *fp, unsigned int matrixSize, unsigned int (*matrix)[matrixSize]);