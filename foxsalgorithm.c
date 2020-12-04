#include "foxsalgorithm.h"
#include <memory.h>
#include <limits.h>
#include <math.h>
#include "matrixutil.h"
#include "parser.h"
#include "matrix_multiplication.h"

#define REORDER 1
#define TAG 37
#define PRINT_RANK 0

int verifyArguments(int processCount, int matrixSize, unsigned int *Q) {

    int maxQ = (int) (sqrt(processCount));

    int possibleProcCount = maxQ * maxQ;

    if (possibleProcCount != processCount) {

        //The number of processes is not a perfect square.

        return 0;
    }

    if (matrixSize % maxQ == 0) {

        *Q = maxQ;

        return 1;
    }

    return 0;
}

int setupGrid(MpiInfo *info) {

    int oldRank;

    int dimensions[2], periods[2], coordinates[2], varying_coords[2];

    //Get the number of processes
    MPI_Comm_size(MPI_COMM_WORLD, &(info->processCount));

    //Get our rank before setting up the grid
    MPI_Comm_rank(MPI_COMM_WORLD, &oldRank);

    info->Q = (int) sqrt(info->processCount);

    //Setup the dimensions of the grid (We want it to be QxQ)
    dimensions[0] = dimensions[1] = info->Q;

    //Specify that the dimensions are circular, not linear
    periods[0] = periods[1] = 1;

    //We want to create a grid with 2 dimensions (Rows and Cols)
    MPI_Cart_create(MPI_COMM_WORLD, 2, dimensions, periods, REORDER, &(info->gridComm));

    //Get our new rank in the grid
#if REORDER
    //We only want to get the new rank if the ranks are changeable (reorder is true)
    MPI_Comm_rank(info->gridComm, &(info->myRank));
#else
    info->myRank = oldRank;
#endif

    //Get our coordinates inside the grid
    MPI_Cart_coords(info->gridComm, info->myRank, 2,
                    coordinates);

    info->myRow = coordinates[0];
    info->myColumn = coordinates[1];

    //Set up row communicator (We want to vary the column, but not the row)
    varying_coords[0] = 0;
    varying_coords[1] = 1;

    MPI_Cart_sub(info->gridComm, varying_coords, &(info->rowComm));

    //Set up column communicator (We want to vary the row, but not the column)
    varying_coords[0] = 1;
    varying_coords[1] = 0;

    MPI_Cart_sub(info->gridComm, varying_coords, &(info->colComm));

    //Successfully setup the grid communicators
    return 1;
}

int setupDatatype(MpiInfo *info, int matrixSize) {

    //The size of the matrix per process
    int perProcessMatrixSize = (matrixSize / info->Q);

    MPI_Type_vector(perProcessMatrixSize, perProcessMatrixSize, perProcessMatrixSize,
                    MPI_UINT32_T, &(info->datatype));

    //Let MPI know about the new datatype
    MPI_Type_commit(&(info->datatype));

    return 1;
}


void performFox(int matrixSize, MpiInfo *info, unsigned int *localA,
                unsigned int *localB, unsigned int *localC) {

    unsigned int n_bar = matrixSize / info->Q;

    int bcast_root,
    //We get the remainder of the division to be able to return to the first process in the row/column if
    //We are the last process in the row/column (Circular shift)
    source = (info->myRow + 1) % info->Q,
            dest = (info->myRow + info->Q - 1) % info->Q;

    unsigned int *tempA = NULL;

    matrix_alloc(n_bar, &tempA);

    matrix_fill(n_bar, tempA, 0);

    for (int step = 0; step < info->Q; step++) {

        //Calculate the chosen column for each row
        bcast_root = (info->myRow + step) % info->Q;

        //If I'm the selected broadcaster from my row
        if (bcast_root == info->myColumn) {

            //Send the localA matrix to the other processes on our row
            MPI_Bcast(localA, 1, info->datatype, bcast_root, info->rowComm);

            doSpecialMatrixMultiply(n_bar, localA, localB, localC);

        } else {
            //Send the tempA matrix to the other processes in our row
            MPI_Bcast(tempA, 1, info->datatype, bcast_root, info->rowComm);

            doSpecialMatrixMultiply(n_bar, tempA, localB, localC);
        }
        //Send our local B matrix to the process directly above

        MPI_Status status;
        //Send our local B matrix to the process directly above
        //Receive the local B matrix from the process directly below
        MPI_Sendrecv_replace(localB, 1, info->datatype, dest, TAG, source, TAG, info->colComm, &status);

//        printf("Status: %d\n", status.MPI_ERROR);

    }

    matrix_free(&tempA);

}


void doAllPairsShortestPathFox(int matrixSize, MpiInfo *info, unsigned int *localA,
                               unsigned int *localB, unsigned int *localC) {
    int m = 1;

    matrix_fill(matrixSize / info->Q, localC, INT_MAX - 1);

    while (m < matrixSize - 1) {
        performFox(matrixSize, info, localA, localB, localC);

        matrix_copy(matrixSize / info->Q, localA, localC);
        matrix_copy(matrixSize / info->Q, localB, localC);

        m *= 2;
    }

}

#define START_ROW(processID, Q, perMatrixSize) ((processID / Q) * perMatrixSize)
#define START_COLUMN(processID, Q, perMatrixSize) ((processID % Q) * perMatrixSize)

int assembleMatrix(unsigned int matrixSize,
                   unsigned int Q,
                   unsigned int processes,
                   unsigned int *originalMatrix,
                   unsigned int *destinationMatrix) {

    unsigned int perMatrix = matrixSize / Q;

    unsigned int k = 0;

    for (int proc = 0; proc < processes; proc++) {

        int row = START_ROW(proc, Q, perMatrix),
                col = START_COLUMN(proc, Q, perMatrix);

        for (int i = 0; i < perMatrix; i++) {
            for (int j = 0; j < perMatrix; j++) {

                destinationMatrix[PROJECT(matrixSize, row + i, col + j)] = originalMatrix[k];
                k++;

            }
        }

    }

    return 1;
}


int buildScatterMatrix(unsigned int matrixSize,
                       unsigned int Q,
                       unsigned int processes,
                       unsigned int *originalMatrix,
                       unsigned int *destinationMatrix) {

    int k = 0;

    unsigned int perMatrixSize = matrixSize / Q;

    for (int proc = 0; proc < processes; proc++) {

        int row = START_ROW(proc, Q, perMatrixSize),
                column = START_COLUMN(proc, Q, perMatrixSize);

        for (int i = 0; i < perMatrixSize; i++) {
            for (int j = 0; j < perMatrixSize; j++) {

                destinationMatrix[k] = originalMatrix[PROJECT(matrixSize, (i + row), (j + column))];

                k++;
            }

        }

    }

    return 1;
}