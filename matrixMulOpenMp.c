#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <immintrin.h>
#include<omp.h>


int **allocateMatrixint(int rows, int cols) {
	int **matrix = (int **)aligned_alloc(32,rows * sizeof(int *));
	for (int i = 0; i < rows; ++i) {
    	matrix[i] = (int *)aligned_alloc(32, cols * sizeof(int));
	}
	return matrix;
}
float **allocateMatrixfloat(int rows, int cols) {
	float **matrix = (float **)aligned_alloc(32,rows * sizeof(float *));
	for (int i = 0; i < rows; ++i) {
    	matrix[i] = (float *)aligned_alloc(32, cols * sizeof(float));
	}
	return matrix;
   
}
long int **allocateMatrixlongint(int rows, int cols) {
	long int **matrix = (long int **)aligned_alloc(32,rows * sizeof(long int *));
	for (int i = 0; i < rows; ++i) {
    	matrix[i] = (long int *)aligned_alloc(32, cols * sizeof(long int));
	}
	return matrix;
   
}
double **allocateMatrixdouble(int rows, int cols) {
	double **matrix = (double **)aligned_alloc(32,rows * sizeof(double *));
	for (int i = 0; i < rows; ++i) {
    	matrix[i] = (double *)aligned_alloc(32, cols * sizeof(double));
	}
	return matrix;
   
}
void readMatrixFromFileint(FILE *file, int **matrix, int rows, int cols) {
	for (int i = 0; i < rows; ++i) {
    	for (int j = 0; j < cols; ++j) {
        	fscanf(file, "%d", &matrix[i][j]);
    	}
	}
}
void readMatrixFromFilefloat(FILE *file, float **matrix, int rows, int cols) {
	for (int i = 0; i < rows; ++i) {
    	for (int j = 0; j < cols; ++j) {
        	fscanf(file, "%f", &matrix[i][j]);
    	}
	}
}
void readMatrixFromFilelongint(FILE *file, long int **matrix, int rows, int cols) {
	for (int i = 0; i < rows; ++i) {
    	for (int j = 0; j < cols; ++j) {
        	fscanf(file, "%ld", &matrix[i][j]);
    	}
	}
}
void readMatrixFromFiledouble(FILE *file, double **matrix, int rows, int cols) {
	for (int i = 0; i < rows; ++i) {
    	for (int j = 0; j < cols; ++j) {
        	fscanf(file, "%lf", &matrix[i][j]);
    	}
	}
}
void multiplyMatricesint(int **matrix_1, int **matrix_2, int **result, int rows1, int cols1, int rows2, int cols2) {
int i,j,k;
#pragma omp parallel for schedule(static) num_threads(4) private(j,k)
	for (int i = 0; i < rows1; i++) {
    	for (int j = 0; j < cols2; j += 8) {
        	__m256i sum = _mm256_setzero_si256();
        	for (int k = 0; k < cols1; k++) {
            	__m256i bc_mat1 = _mm256_set1_epi32(matrix_1[i][k]);
            	__m256i vec_mat2 = _mm256_loadu_si256((__m256i *)&matrix_2[k][j]);
            	__m256i prod = _mm256_mullo_epi32(bc_mat1, vec_mat2);
            	sum = _mm256_add_epi32(sum, prod);
        	}




        	_mm256_storeu_si256((__m256i *)&result[i][j], sum);
    	}
	}
}
void multiplyMatricesfloat(float **matrix_1, float **matrix_2, float **result, int rows1, int cols1, int rows2, int cols2) {
int i,j,k;
#pragma omp parallel for schedule(static) num_threads(4) private(j,k)
for (int i = 0; i < rows1; i++)
	{
    	for (int j = 0; j < cols2; j += 4)
    	{
        	__m256 sum = _mm256_setzero_ps();
        	for (int k = 0; k < cols1; k++)
        	{
            	__m256 bc_mat1 = _mm256_set1_ps(matrix_1[i][k]);
            	__m256 vec_mat2 = _mm256_loadu_ps(&matrix_2[k][j]);
            	__m256 prod = _mm256_mul_ps(bc_mat1, vec_mat2);
            	sum = _mm256_add_ps(sum, prod);
        	}




        	_mm256_storeu_ps(&result[i][j], sum);
    	}
	}
}
void multiplyMatriceslongint(long int **matrix_1, long int **matrix_2, long int **result, int rows1, int cols1, int rows2, int cols2) {
int i,j,k;
#pragma omp parallel for schedule(static) num_threads(4) private(j,k)
for (int i = 0; i < rows1; i++) {
	for (int j = 0; j < cols2; j += 4) {
    	__m256i sum = _mm256_setzero_si256();
    	for (int k = 0; k < cols1; k++) {
        	__m256i bc_mat1 = _mm256_set1_epi64x(matrix_1[i][k]);
        	__m256i vec_mat2 = _mm256_loadu_si256((__m256i *)&matrix_2[k][j]);
        	__m256i prod = _mm256_mul_epi32(bc_mat1, vec_mat2);
        	sum = _mm256_add_epi64(sum, prod);
    	}




    	_mm256_storeu_si256((__m256i *)&result[i][j], sum);  
	}
}

}
void multiplyMatricesdouble(double **matrix_1, double **matrix_2, double **result, int rows1, int cols1, int rows2, int cols2) {
int i,j,k;
#pragma omp parallel for schedule(static) num_threads(4) private(j,k)
for (int i = 0; i < rows1; i++) {
	for (int j = 0; j < cols2; j += 4) {
    	__m256d sum = _mm256_setzero_pd();
    	for (int k = 0; k < cols1; k++) {
        	__m256d bc_mat1 = _mm256_set1_pd(matrix_1[i][k]);
        	__m256d vec_mat2 = _mm256_loadu_pd(&matrix_2[k][j]);
        	__m256d prod = _mm256_mul_pd(bc_mat1, vec_mat2);
        	sum = _mm256_add_pd(sum, prod);
    	}
    	_mm256_storeu_pd(&result[i][j], sum);
	}
}

}

void freeMatrix(int **matrix, int rows) {
	for (int i = 0; i < rows; ++i) {
    	free(matrix[i]);
	}
	free(matrix);
}
void freeMatrixfloat(float **matrix, int rows) {
	for (int i = 0; i < rows; ++i) {
    	free(matrix[i]);
	}
	free(matrix);
}
void freeMatrixlongint(long int **matrix, int rows) {
	for (int i = 0; i < rows; ++i) {
    	free(matrix[i]);
	}
	free(matrix);
}
void freeMatrixdouble(double **matrix, int rows) {
	for (int i = 0; i < rows; ++i) {
    	free(matrix[i]);
	}
	free(matrix);
}
int main(int argc, char *argv[]) {
	clock_t start, end;
	double cpu_time_used;


	FILE *inputFile = fopen(argv[1], "r");
	if (inputFile == NULL) {
    	perror("Error opening input file");
    	return 1;
	}
	int opcode, matrixType;
	int rows1, cols1;
	int rows2, cols2;
	opcode = matrixType = 0;
	fscanf(inputFile, "%d %d", &opcode, &matrixType);
	if (opcode == 4) {//just for checking
    	if (matrixType == 1) {
        	fscanf(inputFile, "%d %d", &rows1, &cols1);
    	 
        	int **matrix1 = allocateMatrixint(rows1, cols1);
        	readMatrixFromFileint(inputFile, matrix1, rows1, cols1);
        	fscanf(inputFile, "%d %d", &rows2, &cols2);
        	int **matrix2 = allocateMatrixint(rows2, cols2);
        	readMatrixFromFileint(inputFile, matrix2, rows2, cols2);
        	fclose(inputFile);
        	if (cols1 != rows2) {
            	printf("Error: Matrices cannot be multiplied.\n");
            	return 1;
        	}
        	int **result = allocateMatrixint(rows1, cols2);
        	start = clock();
        	multiplyMatricesint(matrix1, matrix2, result, rows1, cols1, rows2, cols2);
        	end = clock();
        	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
        	FILE *outputFile = fopen(argv[2], "w");
        	if (outputFile == NULL) {
            	perror("Error opening output file");
            	return 1;
        	}
        	fprintf(outputFile, "%d\n", matrixType);
        	fprintf(outputFile, "%d %d\n", rows1, cols2);
        	for (int i = 0; i < rows1; ++i) {
            	for (int j = 0; j < cols2; ++j) {
                	fprintf(outputFile, "%d ", result[i][j]);
            	}
            	fprintf(outputFile, "\n");
        	}
        	fprintf(outputFile, "Execution Time: %f seconds\n", cpu_time_used);
        	printf("Execution Time: %f seconds\n", cpu_time_used);
        	fclose(outputFile);
        	freeMatrix(matrix1, rows1);
        	freeMatrix(matrix2, rows2);
        	freeMatrix(result, rows1);
    	}
    	else if(matrixType == 2)
    	{
        	fscanf(inputFile, "%d %d", &rows1, &cols1);
        	float **matrix1 = allocateMatrixfloat(rows1, cols1);
        	readMatrixFromFilefloat(inputFile, matrix1, rows1, cols1);
        	fscanf(inputFile, "%d %d", &rows2, &cols2);
        	float **matrix2 = allocateMatrixfloat(rows2, cols2);
        	readMatrixFromFilefloat(inputFile, matrix2, rows2, cols2);
        	fclose(inputFile);
        	if (cols1 != rows2) {
            	printf("Error: Matrices cannot be multiplied.\n");
            	return 1;
        	}
         	float **result = allocateMatrixfloat(rows1, cols2);
        	start = clock();
        	multiplyMatricesfloat(matrix1, matrix2, result, rows1, cols1, rows2, cols2);
        	end = clock();
        	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
        	FILE *outputFile = fopen(argv[2], "w");
        	if (outputFile == NULL) {
            	perror("Error opening output file");
            	return 1;
        	}
        	fprintf(outputFile, "%d\n", matrixType);
        	fprintf(outputFile, "%d %d\n", rows1, cols2);
        	for (int i = 0; i < rows1; ++i) {
            	for (int j = 0; j < cols2; ++j) {
                	fprintf(outputFile, "%f ", result[i][j]);
            	}
            	fprintf(outputFile, "\n");
        	}
        	fprintf(outputFile, "Execution Time: %f seconds\n", cpu_time_used);
        	printf("Execution Time: %f seconds\n", cpu_time_used);
        	fclose(outputFile);
        	freeMatrixfloat(matrix1, rows1);
        	freeMatrixfloat(matrix2, rows2);
        	freeMatrixfloat(result, rows1);
    	}
    	else if(matrixType == 3)
    	{
        	fscanf(inputFile, "%d %d", &rows1, &cols1);
        	long int **matrix1 = allocateMatrixlongint(rows1, cols1);
        	readMatrixFromFilelongint(inputFile, matrix1, rows1, cols1);
        	fscanf(inputFile, "%d %d", &rows2, &cols2);
        	long int **matrix2 = allocateMatrixlongint(rows2, cols2);
        	readMatrixFromFilelongint(inputFile, matrix2, rows2, cols2);
        	fclose(inputFile);
        	if (cols1 != rows2) {
            	printf("Error: Matrices cannot be multiplied.\n");
            	return 1;
        	}
        	long int **result = allocateMatrixlongint(rows1, cols2);
        	start = clock();
        	multiplyMatriceslongint(matrix1, matrix2, result, rows1, cols1, rows2, cols2);
        	end = clock();
        	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;
        	FILE *outputFile = fopen(argv[2], "w");
        	if (outputFile == NULL) {
            	perror("Error opening output file");
            	return 1;
        	}
        	fprintf(outputFile, "%d\n", matrixType);
        	fprintf(outputFile, "%d %d\n", rows1, cols2);
        	for (int i = 0; i < rows1; ++i) {
            	for (int j = 0; j < cols2; ++j) {
                	fprintf(outputFile, "%ld ", result[i][j]);
            	}
            	fprintf(outputFile, "\n");
        	}
        	fprintf(outputFile, "Execution Time: %f seconds\n", cpu_time_used);
        	printf("Execution Time: %f seconds\n", cpu_time_used);
        	fclose(outputFile);
        	freeMatrixlongint(matrix1, rows1);
        	freeMatrixlongint(matrix2, rows2);
        	freeMatrixlongint(result, rows1);
    	}
    	else if(matrixType == 4)
    	{
        	fscanf(inputFile, "%d %d", &rows1, &cols1);
        	double **matrix1 = allocateMatrixdouble(rows1, cols1);
        	readMatrixFromFiledouble(inputFile, matrix1, rows1, cols1);
        	fscanf(inputFile, "%d %d", &rows2, &cols2);
        	double **matrix2 = allocateMatrixdouble(rows2, cols2);
        	readMatrixFromFiledouble(inputFile, matrix2, rows2, cols2);
        	fclose(inputFile);




        	if (cols1 != rows2) {
            	printf("Error: Matrices cannot be multiplied.\n");
            	return 1;
        	}




        	double **result = allocateMatrixdouble(rows1, cols2);
        	start = clock();
        	multiplyMatricesdouble(matrix1, matrix2, result, rows1, cols1, rows2, cols2);
        	end = clock();
        	cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;




        	FILE *outputFile = fopen(argv[2], "w");
        	if (outputFile == NULL) {
            	perror("Error opening output file");
            	return 1;
        	}




        	fprintf(outputFile, "%d\n", matrixType);
        	fprintf(outputFile, "%d %d\n", rows1, cols2);
        	for (int i = 0; i < rows1; ++i) {
            	for (int j = 0; j < cols2; ++j) {
                	fprintf(outputFile, "%lf ", result[i][j]);
            	}
            	fprintf(outputFile, "\n");
        	}




        	fprintf(outputFile, "Execution Time: %f seconds\n", cpu_time_used);
        	printf("Execution Time: %f seconds\n", cpu_time_used);




        	fclose(outputFile);
        	freeMatrixdouble(matrix1, rows1);
        	freeMatrixdouble(matrix2, rows2);
        	freeMatrixdouble(result, rows1);
	 
    	}
    	else
    	{
        	printf("Data Type not supported\n");
    	}
	} else {
    	printf("Opcode not supported\n");
	}




	return 0;
}



































