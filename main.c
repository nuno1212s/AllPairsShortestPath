#include <stdio.h>
#include <stdlib.h>
#include "mpi/mpi.h"
#include "memory.h"
#include <math.h>
#include "parser.h"
#include "foxsalgorithm.h"
#include "matrix_multiplication.h"

int verifyArguments(int processCount, int matrixSize, int *Q);

int main(int argc, char **argv) {

    int matrixSize, numProc = 4, rank = 0, Q;

//    MPI_Init(&argc, &argv);

//    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

//    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        FILE *fp = stdin;

//        printf("Reading matrix...\n");

        fscanf(fp, "%d", &matrixSize);

        if (!verifyArguments(numProc, matrixSize, &Q)) {

//        MPI_Finalize();

            fprintf(stderr,
                    "Failed to initialize the process. The number of processes does not match the size of the matrix.\n");;

            printf("Number of processes: %d . Matrix size: %d .\n", numProc, matrixSize);

            return EXIT_FAILURE;
        }

        //The full matrix, allocated in the stack for faster execution
        //unsigned int matrix[matrixSize][matrixSize];

        //Had to change to this, since the larger inputs ran out of stack memory xD
        unsigned int **matrix = malloc(sizeof(unsigned int) * matrixSize * matrixSize);

        if (parseMatrix(fp, matrixSize, matrix)) {

            printf("Parsed matrix.\n");

            unsigned int **matrixCopy = malloc(sizeof(unsigned int) * matrixSize * matrixSize);

            memcpy(matrixCopy, matrix, sizeof(unsigned int) * matrixSize * matrixSize);

            unsigned int **result = malloc(sizeof(unsigned int) * matrixSize * matrixSize);
            //unsigned int result[matrixSize][matrixSize];

            if (prepareMatrixForAllPairs(matrixSize, matrix)) {

                if (doAllPairsShortestPaths(matrixSize, matrix, result)) {

                    printMatrix(stdout, matrixSize, result);

                }
            }

            free(result);

            //The divided matrixes
//        unsigned int dividedMatrices[numProc][matrixSize / Q][matrixSize / Q];

            unsigned int **dividedMatrix = malloc(sizeof(unsigned int) * (matrixSize / Q) * (matrixSize / Q));

            for (int i = 0; i < numProc; i++) {
                if (subdivideMatrix(matrixSize, matrixCopy, Q, i, dividedMatrix)) {

                    printf("Matrix for process %d:\n", i);

                    printMatrix(stdout, matrixSize / Q, dividedMatrix);
                }

            }

        } else {
            fprintf(stderr, "Failed to parse matrix. \n");
        }

        free(matrix);
    }

    return 0;
}

