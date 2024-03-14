# Optimization with Parallelization.
Read the pdf for answers to questions which were there in the assignment.
# How To run Codes
Just type this into your console to generate a random matrix
```
python Matrixgen.py
```
# Python 
```
python matrix_mul.py
```
# C
```
gcc matrix_mul.c -o cmm
```
# Java
```
javac MatrixMuliplication
java MatrixMultiplication
```
# C with Threads
```
gcc matrixMulThread.c -o mmt -lpthread
```
# C with Vector Instructions
```
gcc -mavx2 -march=haswell matrixMulVector.c -o mv

```
# C with Threads + Vector Instructions
```
gcc -mavx2 -march=haswell matrixMulVectorThread.c -o mvt -lpthread
```
