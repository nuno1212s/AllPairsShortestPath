cmake_minimum_required(VERSION 3.17)
project(Trabalho_1 C)

set(CMAKE_C_STANDARD 11)

add_executable(Trabalho_1 main.c parser.c parser.h foxsalgorithm.c foxsalgorithm.h matrix_multiplication.c matrix_multiplication.h)
target_link_libraries(Trabalho_1 m)