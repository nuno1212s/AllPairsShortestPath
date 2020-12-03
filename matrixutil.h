#ifndef TRABALHO_1_MATRIXUTIL_H
#define TRABALHO_1_MATRIXUTIL_H

#define PROJECT(size, i, j) (((i) * size) + j)

int matrix_alloc(unsigned int size, unsigned int **matrix);

int matrix_copy(unsigned int size, unsigned int *dest, unsigned int *source);

int matrix_fill(unsigned int matrixSize, unsigned int *matrix, unsigned int value);

int matrix_free(unsigned int **matrix);

#endif //TRABALHO_1_MATRIXUTIL_H
