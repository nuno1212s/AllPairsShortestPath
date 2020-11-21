
#include <mpi/mpi.h>

struct MpiInfo_ {
    int processCount; //Total amount of processes

    MPI_Comm gridComm; //The communicator for the entire grid (All processes)
    MPI_Comm rowComm; //The communicator for the row of this process
    MPI_Comm colComm; //The communicator for the column of this process

    int Q; //The Q value (Q * Q = ProcessCount)

    int myRow, myColumn, myRank; //The row, column and rank of this process

    MPI_Datatype datatype;
};

typedef struct MpiInfo_ MpiInfo;

/**
 * Set up the communicators, read the matrix and send the input to all processes
 *
 * @param info The allocated info, to be filled out by this method
 * @return
 */
int setupGrid(MpiInfo *info);

/**
 * Setup the datatype needed for MPI
 *
 * This will place the Datatype into the MpiInfo struct
 *
 * @param info The MPI Info
 * @param matrixSize The size of the matrix we are going to process
 * @return
 */
int setupDatatype(MpiInfo *info, int matrixSize);

/**
 * Check if the arguments are correct (Number of processes correct for the size of the matrix that we have)
 *
 * @param processCount The amount of processes
 * @param matrixSize The size of the matrix
 * @param Q
 * @return
 */
int verifyArguments(int processCount, int matrixSize, int *Q);

void performFox(int matrixSize, MpiInfo *info, unsigned int(*localA)[matrixSize / info->Q],
                unsigned int (*localB)[matrixSize / info->Q], unsigned int (*localC)[matrixSize / info->Q]);

/**
 * Subdivide the matrix and store the matrix for the process processID in the argument destinationMatrix
 *
 * @param matrixSize The size of the original matrix
 * @param matrix The original matrix
 * @param Q The Q of the problem (Sqrt of the number of processes)
 * @param processID The ID of the process we want to get the matrix for
 * @param destinationMatrix Where to put the processes matrix
 * @return
 */
int subdivideMatrix(unsigned int matrixSize, unsigned int (*matrix)[matrixSize],
                    unsigned int Q,
                    unsigned int processID,
                    unsigned int (*destinationMatrix)[matrixSize / Q]);