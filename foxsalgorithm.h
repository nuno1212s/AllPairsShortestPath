
int subdivideMatrix(unsigned int matrixSize, unsigned int (*matrix)[matrixSize],
                unsigned int processCount, unsigned int Q,
                unsigned int (*matrixDestinations)[matrixSize / Q][matrixSize / Q]);