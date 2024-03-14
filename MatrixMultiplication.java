import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintWriter;
import java.util.Scanner;

public class MatrixMultiplication {

    // Function to allocate memory for a matrix
    static int[][] allocateMatrixInt(int rows, int cols) {
        return new int[rows][cols];
    }

    static float[][] allocateMatrixFloat(int rows, int cols) {
        return new float[rows][cols];
    }

    static long[][] allocateMatrixLong(int rows, int cols) {
        return new long[rows][cols];
    }

    static double[][] allocateMatrixDouble(int rows, int cols) {
        return new double[rows][cols];
    }

    // Function to read matrix data from file
    static void readMatrixFromFileInt(Scanner scanner, int[][] matrix, int rows, int cols) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] = scanner.nextInt();
            }
        }
    }

    static void readMatrixFromFileFloat(Scanner scanner, float[][] matrix, int rows, int cols) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] = scanner.nextFloat();
            }
        }
    }

    static void readMatrixFromFileLong(Scanner scanner, long[][] matrix, int rows, int cols) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] = scanner.nextLong();
            }
        }
    }

    static void readMatrixFromFileDouble(Scanner scanner, double[][] matrix, int rows, int cols) {
        for (int i = 0; i < rows; ++i) {
            for (int j = 0; j < cols; ++j) {
                matrix[i][j] = scanner.nextDouble();
            }
        }
    }

    // Function to multiply two matrices
    static void multiplyMatricesInt(int[][] matrix1, int[][] matrix2, int[][] result, int rows1, int cols1, int rows2, int cols2) {
        for (int i = 0; i < rows1; ++i) {
            for (int j = 0; j < cols2; ++j) {
                result[i][j] = 0;
                for (int k = 0; k < cols1; ++k) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
    }

    static void multiplyMatricesFloat(float[][] matrix1, float[][] matrix2, float[][] result, int rows1, int cols1, int rows2, int cols2) {
        for (int i = 0; i < rows1; ++i) {
            for (int j = 0; j < cols2; ++j) {
                result[i][j] = 0;
                for (int k = 0; k < cols1; ++k) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
    }

    static void multiplyMatricesLong(long[][] matrix1, long[][] matrix2, long[][] result, int rows1, int cols1, int rows2, int cols2) {
        for (int i = 0; i < rows1; ++i) {
            for (int j = 0; j < cols2; ++j) {
                result[i][j] = 0;
                for (int k = 0; k < cols1; ++k) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
    }

    static void multiplyMatricesDouble(double[][] matrix1, double[][] matrix2, double[][] result, int rows1, int cols1, int rows2, int cols2) {
        for (int i = 0; i < rows1; ++i) {
            for (int j = 0; j < cols2; ++j) {
                result[i][j] = 0;
                for (int k = 0; k < cols1; ++k) {
                    result[i][j] += matrix1[i][k] * matrix2[k][j];
                }
            }
        }
    }

    // Function to free memory allocated for a matrix
    static void freeMatrixInt(int[][] matrix) {
        // Java garbage collector will take care of this
    }

    static void freeMatrixFloat(float[][] matrix) {
        // Java garbage collector will take care of this
    }

    static void freeMatrixLong(long[][] matrix) {
        // Java garbage collector will take care of this
    }

    static void freeMatrixDouble(double[][] matrix) {
        // Java garbage collector will take care of this
    }

    public static void main(String[] args) {
        long start, end;
        double cpuTimeUsed;

        // Open the input file
        try (Scanner inputFile = new Scanner(new File(args[0]))) {
            // Read the opcode and matrix type
            int opcode = 0, matrixType = 0;
            int rows1, cols1;
            int rows2, cols2;

            opcode = matrixType = 0;
            opcode = inputFile.nextInt();

            if (opcode == 3) {
                matrixType = inputFile.nextInt();

                if (matrixType == 1) {
                    rows1 = inputFile.nextInt();
                    cols1 = inputFile.nextInt();
                    int[][] matrix1 = allocateMatrixInt(rows1, cols1);
                    readMatrixFromFileInt(inputFile, matrix1, rows1, cols1);

                    rows2 = inputFile.nextInt();
                    cols2 = inputFile.nextInt();
                    int[][] matrix2 = allocateMatrixInt(rows2, cols2);
                    readMatrixFromFileInt(inputFile, matrix2, rows2, cols2);

                    inputFile.close();

                    if (cols1 != rows2) {
                        System.out.println("Error: Matrices cannot be multiplied.");
                        return;
                    }

                    int[][] result = allocateMatrixInt(rows1, cols2);

                    start = System.currentTimeMillis();
                    multiplyMatricesInt(matrix1, matrix2, result, rows1, cols1, rows2, cols2);
                    end = System.currentTimeMillis();
                    cpuTimeUsed = (double) (end - start) / 1000;

                    try (PrintWriter outputFile = new PrintWriter(new File(args[1]))) {
                        outputFile.println(matrixType);
                        outputFile.println(rows1 + " " + cols2);

                        for (int i = 0; i < rows1; ++i) {
                            for (int j = 0; j < cols2; ++j) {
                                outputFile.print(result[i][j] + " ");
                            }
                            outputFile.println();
                        }

                        outputFile.println("Execution Time: " + cpuTimeUsed + " seconds");
                        System.out.println("Execution Time: " + cpuTimeUsed + " seconds");
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                    freeMatrixInt(matrix1);
                    freeMatrixInt(matrix2);
                    freeMatrixInt(result);
                } else if (matrixType == 2) {
                    rows1 = inputFile.nextInt();
                    cols1 = inputFile.nextInt();
                    float[][] matrix1 = allocateMatrixFloat(rows1, cols1);
                    readMatrixFromFileFloat(inputFile, matrix1, rows1, cols1);

                    rows2 = inputFile.nextInt();
                    cols2 = inputFile.nextInt();
                    float[][] matrix2 = allocateMatrixFloat(rows2, cols2);
                    readMatrixFromFileFloat(inputFile, matrix2, rows2, cols2);

                    inputFile.close();

                    if (cols1 != rows2) {
                        System.out.println("Error: Matrices cannot be multiplied.");
                        return;
                    }

                    float[][] result = allocateMatrixFloat(rows1, cols2);

                    start = System.currentTimeMillis();
                    multiplyMatricesFloat(matrix1, matrix2, result, rows1, cols1, rows2, cols2);
                    end = System.currentTimeMillis();
                    cpuTimeUsed = (double) (end - start) / 1000;

                    try (PrintWriter outputFile = new PrintWriter(new File(args[1]))) {
                        outputFile.println(matrixType);
                        outputFile.println(rows1 + " " + cols2);

                        for (int i = 0; i < rows1; ++i) {
                            for (int j = 0; j < cols2; ++j) {
                                outputFile.print(result[i][j] + " ");
                            }
                            outputFile.println();
                        }

                        outputFile.println("Execution Time: " + cpuTimeUsed + " seconds");
                        System.out.println("Execution Time: " + cpuTimeUsed + " seconds");
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                    freeMatrixFloat(matrix1);
                    freeMatrixFloat(matrix2);
                    freeMatrixFloat(result);
                } else if (matrixType == 3) {
                    rows1 = inputFile.nextInt();
                    cols1 = inputFile.nextInt();
                    long[][] matrix1 = allocateMatrixLong(rows1, cols1);
                    readMatrixFromFileLong(inputFile, matrix1, rows1, cols1);

                    rows2 = inputFile.nextInt();
                    cols2 = inputFile.nextInt();
                    long[][] matrix2 = allocateMatrixLong(rows2, cols2);
                    readMatrixFromFileLong(inputFile, matrix2, rows2, cols2);

                    inputFile.close();

                    if (cols1 != rows2) {
                        System.out.println("Error: Matrices cannot be multiplied.");
                        return;
                    }

                    long[][] result = allocateMatrixLong(rows1, cols2);

                    start = System.currentTimeMillis();
                    multiplyMatricesLong(matrix1, matrix2, result, rows1, cols1, rows2, cols2);
                    end = System.currentTimeMillis();
                    cpuTimeUsed = (double) (end - start) / 1000;

                    try (PrintWriter outputFile = new PrintWriter(new File(args[1]))) {
                        outputFile.println(matrixType);
                        outputFile.println(rows1 + " " + cols2);

                        for (int i = 0; i < rows1; ++i) {
                            for (int j = 0; j < cols2; ++j) {
                                outputFile.print(result[i][j] + " ");
                            }
                            outputFile.println();
                        }

                        outputFile.println("Execution Time: " + cpuTimeUsed + " seconds");
                        System.out.println("Execution Time: " + cpuTimeUsed + " seconds");
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                    freeMatrixLong(matrix1);
                    freeMatrixLong(matrix2);
                    freeMatrixLong(result);
                } else if (matrixType == 4) {
                    rows1 = inputFile.nextInt();
                    cols1 = inputFile.nextInt();
                    double[][] matrix1 = allocateMatrixDouble(rows1, cols1);
                    readMatrixFromFileDouble(inputFile, matrix1, rows1, cols1);

                    rows2 = inputFile.nextInt();
                    cols2 = inputFile.nextInt();
                    double[][] matrix2 = allocateMatrixDouble(rows2, cols2);
                    readMatrixFromFileDouble(inputFile, matrix2, rows2, cols2);

                    inputFile.close();

                    if (cols1 != rows2) {
                        System.out.println("Error: Matrices cannot be multiplied.");
                        return;
                    }

                    double[][] result = allocateMatrixDouble(rows1, cols2);

                    start = System.currentTimeMillis();
                    multiplyMatricesDouble(matrix1, matrix2, result, rows1, cols1, rows2, cols2);
                    end = System.currentTimeMillis();
                    cpuTimeUsed = (double) (end - start) / 1000;

                    try (PrintWriter outputFile = new PrintWriter(new File(args[1]))) {
                        outputFile.println(matrixType);
                        outputFile.println(rows1 + " " + cols2);

                        for (int i = 0; i < rows1; ++i) {
                            for (int j = 0; j < cols2; ++j) {
                                outputFile.print(result[i][j] + " ");
                            }
                            outputFile.println();
                        }

                        outputFile.println("Execution Time: " + cpuTimeUsed + " seconds");
                        System.out.println("Execution Time: " + cpuTimeUsed + " seconds");
                    } catch (FileNotFoundException e) {
                        e.printStackTrace();
                    }

                    freeMatrixDouble(matrix1);
                    freeMatrixDouble(matrix2);
                    freeMatrixDouble(result);
                } else {
                    System.out.println("Data Type not supported");
                }
            } else {
                System.out.println("Opcode not supported");
            }

        } catch (FileNotFoundException e) {
            e.printStackTrace();
        }
    }
}
