#include <stdio.h>
#include <stdlib.h>
#include <zconf.h>
#include "mpi/mpi.h"
#include "parser.h"
#include "foxsalgorithm.h"
#include "matrix_multiplication.h"
#include "matrixutil.h"

static void resetMatrix(unsigned int matrixSize, unsigned int (*toReset)[matrixSize]) {

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            toReset[i][j] = 0;
        }
    }

}

int main(int argc, char **argv) {


    int matrixSize, numProc = 4, rank = 0, Q;

    MPI_Init(&argc, &argv);

    MPI_Comm_size(MPI_COMM_WORLD, &numProc);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    FILE *fp = stdin;

    if (rank == 0) {
        fscanf(fp, "%d", &matrixSize);
    }

    MPI_Bcast(&matrixSize, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    unsigned int *dividedMatrix = NULL;

    if (rank == 0) {

        matrix_alloc(matrixSize, &dividedMatrix);

        if (!verifyArguments(numProc, matrixSize, &Q)) {

            MPI_Finalize();

            fprintf(stderr,
                    "Failed to initialize the process. The number of processes does not match the size of the matrix.\n");;

            printf("Number of processes: %d . Matrix size: %d .\n", numProc, matrixSize);

            return EXIT_FAILURE;
        }

        //The full matrix, allocated in the stack for faster execution
        //unsigned int matrix[matrixSize][matrixSize];

        //Had to change to this, since the larger inputs ran out of stack memory xD
        unsigned int *matrix = NULL;

        matrix_alloc(matrixSize, &matrix);

        if (parseMatrix(fp, matrixSize, matrix)) {

            prepareMatrixForAllPairs(matrixSize, matrix);

            printf("Parsed matrix.\n");

            unsigned int perProcessMatrixSize = (matrixSize / Q);

            printf("Initial matrix:\n");

            printMatrix(stdout, matrixSize, matrix);

            buildScatterMatrix(matrixSize, Q, numProc, matrix, dividedMatrix);

            printf("Scatter matrix: \n");

            printMatrix(stdout, matrixSize, dividedMatrix);


        } else {
            fprintf(stderr, "Failed to parse matrix. \n");

            exit(EXIT_FAILURE);
        }

        matrix_free(&matrix);
    }

//    if (rank == 0) {
//        volatile int i = 0;
//        char hostname[256];
//        gethostname(hostname, sizeof(hostname));
//        printf("PID %d on %s ready for attach\n", getpid(), hostname);
//        fflush(stdout);
//        while (0 == i)
//            sleep(5);
//    }

    MPI_Bcast(&Q, 1, MPI_UINT32_T, 0, MPI_COMM_WORLD);

    MpiInfo info;

    setupGrid(&info);

    setupDatatype(&info, matrixSize);

//    printf("INFO: Q: %d, Process Count: %d, My Row: %d, My Column: %d, My Rank: %d \n", info.Q,
//           info.processCount, info.myRow, info.myColumn, info.myRank);

    unsigned int *localA = NULL, *localB = NULL, *localC = NULL;

    matrix_alloc(matrixSize / Q, &localA);
    matrix_alloc(matrixSize / Q, &localB);
    matrix_alloc(matrixSize / Q, &localC);

    MPI_Scatter(dividedMatrix, 1, info.datatype, localA, 1, info.datatype, 0, info.gridComm);

    matrix_copy(matrixSize / Q, localB, localA);

    if (info.myRank == 0) {
        printf("Process %d received: \n Local A:\n", rank);

        printMatrix(stdout, matrixSize / Q, localA);

        printf("LocalB: \n");

        printMatrix(stdout, matrixSize / Q, localB);
    }

    doAllPairsShortestPathFox(matrixSize, &info, localA, localB, localC);

    if (info.myRank == 0) {
        printf("Process %d ended with: \n", rank);

        printMatrix(stdout, matrixSize / Q, localC);
        printf("\n");
    }

    MPI_Gather(localC, 1, info.datatype, dividedMatrix, 1, info.datatype, 0, info.gridComm);

    if (info.myRank == 0) {
        printf("Final matrix:\n");

        printMatrix(stdout, matrixSize, dividedMatrix);
    }

    if (dividedMatrix != NULL)
        free(dividedMatrix);

    free(localA);
    free(localB);
    free(localC);

    return 0;
}

