cmake_minimum_required(VERSION 3.17)
project(Trabalho_1 C)
find_package(MPI REQUIRED)

set(CMAKE_C_STANDARD 11)
SET(CMAKE_C_COMPILER "mpicc")
SET(CMAKE_CXX_COMPILER "mpic++")

add_executable(Trabalho_1 main.c parser.c parser.h foxsalgorithm.c foxsalgorithm.h matrix_multiplication.c matrix_multiplication.h matrixutil.c matrixutil.h)
target_link_libraries(Trabalho_1 m)