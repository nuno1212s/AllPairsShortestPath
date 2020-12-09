#include <stdio.h>
#include <stdlib.h>
#include "mpi.h"
#include "parser.h"
#include "foxsalgorithm.h"
#include "matrix_multiplication.h"
#include "matrixutil.h"

#define ROOT 0

int main(int argc, char **argv) {

    int numProc = 4, rank = 0;

    unsigned int matrixSize;
    unsigned int Q;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FILE *fp = stdin;

    if (rank == 0) {
        fscanf(fp, "%d", &matrixSize);
    }

    MPI_Bcast(&matrixSize, 1, MPI_UINT32_T, ROOT, MPI_COMM_WORLD);

    unsigned int *dividedMatrix = NULL;

    matrix_alloc(matrixSize, &dividedMatrix);

    if (rank == 0) {

        if (!verifyArguments(numProc, matrixSize, &Q)) {

            fprintf(stderr,
                    "Failed to initialize the process. The number of processes is not a perfect square or does not match the size of the matrix.\n");;

            printf("Number of processes: %d . Matrix size: %d .\n", numProc, matrixSize);

            matrix_free(&dividedMatrix);

            exit(EXIT_FAILURE);

            return EXIT_FAILURE;
        }

        //The full matrix, allocated in the stack for faster execution
        //unsigned int matrix[matrixSize][matrixSize];

        //Had to change to this, since the larger inputs ran out of stack memory xD
        unsigned int *matrix = NULL;

        matrix_alloc(matrixSize, &matrix);

        if (parseMatrix(fp, matrixSize, matrix)) {

            prepareMatrixForAllPairs(matrixSize, matrix);

            buildScatterMatrix(matrixSize, Q, numProc, matrix, dividedMatrix);

        } else {
            fprintf(stderr, "Failed to parse matrix. \n");

            exit(EXIT_FAILURE);
        }

        matrix_free(&matrix);
    }

    MPI_Bcast(&Q, 1, MPI_UINT32_T, ROOT, MPI_COMM_WORLD);

    MpiInfo info;

    setupGrid(&info);

    setupDatatype(&info, matrixSize);

    unsigned int *localA = NULL, *localB = NULL, *localC = NULL;

    matrix_alloc(matrixSize / Q, &localA);
    matrix_alloc(matrixSize / Q, &localB);
    matrix_alloc(matrixSize / Q, &localC);

    MPI_Scatter(dividedMatrix, 1, info.datatype, localA, 1, info.datatype, ROOT, info.gridComm);

    matrix_copy(matrixSize / Q, localB, localA);

    doAllPairsShortestPathFox(matrixSize, &info, localA, localB, localC);

    MPI_Gather(localC, 1, info.datatype, dividedMatrix, 1, info.datatype, ROOT, info.gridComm);

    if (info.myRank == ROOT) {
        unsigned int *finalDestination = NULL;

        matrix_alloc(matrixSize, &finalDestination);

        assembleMatrix(matrixSize, Q, info.processCount, dividedMatrix, finalDestination);

        FILE *resultFile = stdout;

        printMatrix(resultFile, matrixSize, finalDestination);

        matrix_free(&finalDestination);
    }

    matrix_free(&dividedMatrix);

    matrix_free(&localA);
    matrix_free(&localB);
    matrix_free(&localC);

    MPI_Type_free(&info.datatype);

    MPI_Finalize();

    return 0;
}

