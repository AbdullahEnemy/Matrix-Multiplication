#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <time.h>
#include <immintrin.h>


#define MAX_THREADS 8


typedef struct {
	int thread_id;
	int start_row;
	int end_row;
	int **matrix1;
	int **matrix2;
	int **result;
	int rows1;
	int cols1;
	int rows2;
	int cols2;
} ThreadArgsint;
typedef struct {
	int thread_id;
	int start_row;
	int end_row;
	float **matrix1;
	float **matrix2;
	float **result;
	int rows1;
	int cols1;
	int rows2;
	int cols2;
} ThreadArgsfloat;
typedef struct {
	int thread_id;
	int start_row;
	int end_row;
	long int **matrix1;
	long int **matrix2;
	long int **result;
	int rows1;
	int cols1;
	int rows2;
	int cols2;
} ThreadArgslongint;
typedef struct {
	int thread_id;
	int start_row;
	int end_row;
	double **matrix1;
	double **matrix2;
	double **result;
	int rows1;
	int cols1;
	int rows2;
	int cols2;
} ThreadArgsdouble;


int **allocateMatrix(int rows, int cols) {
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


void readMatrixFromFile(FILE *file, int **matrix, int rows, int cols) {
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




void *multiplyMatricesint(void *args) {
	ThreadArgsint *thread_args = (ThreadArgsint *)args;


	for (int i = thread_args->start_row; i < thread_args->end_row; ++i) {
    	for (int j = 0; j < thread_args->cols2; j += 8) {  // Assuming 8 is the vector width (AVX2)
        	__m256i sum = _mm256_setzero_si256();
        	for (int k = 0; k < thread_args->cols1; ++k) {
            	__m256i a = _mm256_set1_epi32(thread_args->matrix1[i][k]);
            	__m256i b = _mm256_loadu_si256((__m256i*)&thread_args->matrix2[k][j]);
            	sum = _mm256_add_epi32(sum, _mm256_mullo_epi32(a, b));
        	}
        	_mm256_storeu_si256((__m256i*)&thread_args->result[i][j], sum);
    	}
	}


	pthread_exit(NULL);
}


void *multiplyMatricesfloat(void *args) {
	ThreadArgsfloat *thread_args = (ThreadArgsfloat *)args;


	for (int i = thread_args->start_row; i < thread_args->end_row; ++i) {
    	for (int j = 0; j < thread_args->cols2; j += 8) {  // Assuming 8 is the vector width (AVX2)
        	__m256 sum = _mm256_setzero_ps();
        	for (int k = 0; k < thread_args->cols1; ++k) {
            	__m256 a = _mm256_set1_ps(thread_args->matrix1[i][k]);
            	__m256 b = _mm256_loadu_ps(&thread_args->matrix2[k][j]);
            	sum = _mm256_add_ps(sum, _mm256_mul_ps(a, b));
        	}
        	_mm256_storeu_ps(&thread_args->result[i][j], sum);
    	}
	}


	pthread_exit(NULL);
}
void *multiplyMatriceslongint(void *args) {
	ThreadArgslongint *thread_args = (ThreadArgslongint *)args;




for (int i = thread_args->start_row; i < thread_args->end_row; i++) {
	for (int j = 0; j < thread_args->cols2; j += 4) {
    	__m256i sum = _mm256_setzero_si256();
    	for (int k = 0; k < thread_args->cols1; k++) {
        	__m256i bc_mat1 = _mm256_set1_epi64x(thread_args->matrix1[i][k]);
        	__m256i vec_mat2 = _mm256_loadu_si256((__m256i*)&thread_args->matrix2[k][j]);
        	__m256i prod = _mm256_mul_epi32(bc_mat1, vec_mat2);
        	sum = _mm256_add_epi64(sum, prod);
    	}




    	_mm256_storeu_si256((__m256i*)&thread_args->result[i][j], sum);  
	}
}




	pthread_exit(NULL);
}












void *multiplyMatricesdouble(void *args) {
	ThreadArgsdouble *thread_args = (ThreadArgsdouble *)args;


	for (int i = thread_args->start_row; i < thread_args->end_row; ++i) {
    	for (int j = 0; j < thread_args->cols2; j += 4) {  // Assuming 4 is the vector width (AVX2)
        	__m256d sum = _mm256_setzero_pd();
        	for (int k = 0; k < thread_args->cols1; ++k) {
            	__m256d a = _mm256_set1_pd(thread_args->matrix1[i][k]);
            	__m256d b = _mm256_loadu_pd(&thread_args->matrix2[k][j]);
            	sum = _mm256_add_pd(sum, _mm256_mul_pd(a, b));
        	}
        	_mm256_storeu_pd(&thread_args->result[i][j], sum);
    	}
	}


	pthread_exit(NULL);
}






int main(int argc, char *argv[]) {
	struct timespec start, end;
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


	if (opcode == 6) {
    	if (matrixType == 1) {
        	fscanf(inputFile, "%d %d", &rows1, &cols1);
        	int **matrix1 = allocateMatrix(rows1, cols1);
        	readMatrixFromFile(inputFile, matrix1, rows1, cols1);
        	fscanf(inputFile, "%d %d", &rows2, &cols2);
        	int **matrix2 = allocateMatrix(rows2, cols2);
        	readMatrixFromFile(inputFile, matrix2, rows2, cols2);
        	fclose(inputFile);


        	if (cols1 != rows2) {
            	printf("Error: Matrices cannot be multiplied.\n");
            	return 1;
        	}


        	int **result = allocateMatrix(rows1, cols2);


        	pthread_t threads[MAX_THREADS];
        	ThreadArgsint thread_args[MAX_THREADS];


        	int rows_per_thread = rows1 / MAX_THREADS;


        	clock_gettime(CLOCK_MONOTONIC, &start);


        	for (int i = 0; i < MAX_THREADS; ++i) {
            	thread_args[i].thread_id = i;
            	thread_args[i].start_row = i * rows_per_thread;
            	thread_args[i].end_row = (i == MAX_THREADS - 1) ? rows1 : (i + 1) * rows_per_thread;
            	thread_args[i].matrix1 = matrix1;
            	thread_args[i].matrix2 = matrix2;
            	thread_args[i].result = result;
            	thread_args[i].rows1 = rows1;
            	thread_args[i].cols1 = cols1;
            	thread_args[i].rows2 = rows2;
            	thread_args[i].cols2 = cols2;


            	pthread_create(&threads[i], NULL, multiplyMatricesint, (void *)&thread_args[i]);
        	}


        	for (int i = 0; i < MAX_THREADS; ++i) {
            	pthread_join(threads[i], NULL);
        	}


        	clock_gettime(CLOCK_MONOTONIC, &end);
        	cpu_time_used = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;


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


        	pthread_t threads[MAX_THREADS];
        	ThreadArgsfloat thread_args[MAX_THREADS];


        	int rows_per_thread = rows1 / MAX_THREADS;


        	clock_gettime(CLOCK_MONOTONIC, &start);


        	for (int i = 0; i < MAX_THREADS; ++i) {
            	thread_args[i].thread_id = i;
            	thread_args[i].start_row = i * rows_per_thread;
            	thread_args[i].end_row = (i == MAX_THREADS - 1) ? rows1 : (i + 1) * rows_per_thread;
            	thread_args[i].matrix1 = matrix1;
            	thread_args[i].matrix2 = matrix2;
            	thread_args[i].result = result;
            	thread_args[i].rows1 = rows1;
            	thread_args[i].cols1 = cols1;
            	thread_args[i].rows2 = rows2;
            	thread_args[i].cols2 = cols2;


            	pthread_create(&threads[i], NULL, multiplyMatricesfloat, (void *)&thread_args[i]);
        	}


        	for (int i = 0; i < MAX_THREADS; ++i) {
            	pthread_join(threads[i], NULL);
        	}


        	clock_gettime(CLOCK_MONOTONIC, &end);
        	cpu_time_used = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;


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




        	pthread_t threads[MAX_THREADS];
        	ThreadArgslongint thread_args[MAX_THREADS];




        	int rows_per_thread = rows1 / MAX_THREADS;




        	clock_gettime(CLOCK_MONOTONIC, &start);




        	for (int i = 0; i < MAX_THREADS; ++i) {
            	thread_args[i].thread_id = i;
            	thread_args[i].start_row = i * rows_per_thread;
            	thread_args[i].end_row = (i == MAX_THREADS - 1) ? rows1 : (i + 1) * rows_per_thread;
            	thread_args[i].matrix1 = matrix1;
            	thread_args[i].matrix2 = matrix2;
            	thread_args[i].result = result;
            	thread_args[i].rows1 = rows1;
            	thread_args[i].cols1 = cols1;
            	thread_args[i].rows2 = rows2;
            	thread_args[i].cols2 = cols2;




            	pthread_create(&threads[i], NULL, multiplyMatriceslongint, (void *)&thread_args[i]);
        	}




        	for (int i = 0; i < MAX_THREADS; ++i) {
            	pthread_join(threads[i], NULL);
        	}




        	clock_gettime(CLOCK_MONOTONIC, &end);
        	cpu_time_used = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;




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




        	pthread_t threads[MAX_THREADS];
        	ThreadArgsdouble thread_args[MAX_THREADS];




        	int rows_per_thread = rows1 / MAX_THREADS;




        	clock_gettime(CLOCK_MONOTONIC, &start);




        	for (int i = 0; i < MAX_THREADS; ++i) {
            	thread_args[i].thread_id = i;
            	thread_args[i].start_row = i * rows_per_thread;
            	thread_args[i].end_row = (i == MAX_THREADS - 1) ? rows1 : (i + 1) * rows_per_thread;
            	thread_args[i].matrix1 = matrix1;
            	thread_args[i].matrix2 = matrix2;
            	thread_args[i].result = result;
            	thread_args[i].rows1 = rows1;
            	thread_args[i].cols1 = cols1;
            	thread_args[i].rows2 = rows2;
            	thread_args[i].cols2 = cols2;




            	pthread_create(&threads[i], NULL, multiplyMatricesdouble, (void *)&thread_args[i]);
        	}

        	for (int i = 0; i < MAX_THREADS; ++i) {
            	pthread_join(threads[i], NULL);
        	}

        	clock_gettime(CLOCK_MONOTONIC, &end);
        	cpu_time_used = (end.tv_sec - start.tv_sec) + (end.tv_nsec - start.tv_nsec) / 1e9;
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
    	printf("Matrixx Type not supported\n");
    	}
	} else {
    	printf("Opcode not supported\n");
	}


	return 0;
}











