#include "foxsalgorithm.h"
#include <stdlib.h>
#include "parser.h"
#include "math.h"
#include "matrix_multiplication.h"

#define REORDER 1
#define TAG 37

static void resetMatrix(unsigned int matrixSize, unsigned int (*toReset)[matrixSize]) {

    for (int i = 0; i < matrixSize; i++) {
        for (int j = 0; j < matrixSize; j++) {
            toReset[i][j] = 0;
        }
    }

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
                    MPI_INT32_T, &(info->datatype));

    //Let MPI know about the new datatype
    MPI_Type_commit(&(info->datatype));

}


void performFox(int matrixSize, MpiInfo *info, unsigned int(*localA)[matrixSize / info->Q],
                unsigned int (*localB)[matrixSize / info->Q], unsigned int (*localC)[matrixSize / info->Q]) {

    unsigned int n_bar = matrixSize / info->Q;

    int bcast_root,
            //We get the remainder of the division to be able to return to the first process in the row/column if
            //We are the last process in the row/column (Circular shift)
            source = (info->myRow + 1) % info->Q,
            dest = (info->myRow + info->Q - 1) % info->Q;

    MPI_Status status;

    resetMatrix(n_bar, localC);

    unsigned int (*tempA)[n_bar] = malloc(sizeof(unsigned int) * n_bar * n_bar);

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
        MPI_Send(localB, 1, info->datatype, dest, TAG, info->colComm);
        //Receive the local B matrix from the process directly below
        MPI_Recv(localB, 1, info->datatype, source, TAG, info->colComm, &status);

    }

    free(tempA);

}

int subdivideMatrix(unsigned int matrixSize, unsigned int (*matrix)[matrixSize],
                    unsigned int Q,
                    unsigned int processID,
                    unsigned int (*matrixDestinations)[matrixSize / Q]) {

    unsigned int perMatrixSize = matrixSize / Q;

    //Omega (x) = (x / Q, x mod Q)
    unsigned int startRow = (processID / Q) * perMatrixSize, endRow = startRow + perMatrixSize,
            startColumn = (processID % Q) * perMatrixSize, endColumn = startColumn + perMatrixSize;

    for (unsigned int i = startRow; i < endRow; i++) {
        for (unsigned int j = startColumn; j < endColumn; j++) {

            matrixDestinations[i - startRow][j - startColumn] = matrix[i][j];

        }
    }

    return 1;
}