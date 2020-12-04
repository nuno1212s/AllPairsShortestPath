CC=mpicc
ARGS=-Wall
LINKS=-lm
OUTPUT=T1

all:
	$(CC) $(ARGS) main.c matrix_multiplication.c matrixutil.c parser.c foxsalgorithm.c -o $(OUTPUT) $(LINKS)

clean:
	rm -f *.o $(OUTPUT)