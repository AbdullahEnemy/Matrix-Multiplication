import time

def read_input_file(filename):
    with open(filename, 'r') as file:
        opcode = int(file.readline().strip())
        data_type = int(file.readline().strip())
        first_matrix_dims = tuple(map(int, file.readline().strip().split()))
        first_matrix = []
        for _ in range(first_matrix_dims[0]):
            first_matrix.append(list(map(lambda x: int(x) if data_type <= 2 else float(x), file.readline().strip().split())))
        second_matrix_dims = tuple(map(int, file.readline().strip().split()))
        second_matrix = []
        for _ in range(second_matrix_dims[0]):
            second_matrix.append(list(map(lambda x: int(x) if data_type <= 2 else float(x), file.readline().strip().split())))
    return opcode, data_type, first_matrix, second_matrix

def multiply_matrices(first_matrix, second_matrix):
    result = [[0 for _ in range(len(second_matrix[0]))] for _ in range(len(first_matrix))]
    for i in range(len(first_matrix)):
        for j in range(len(second_matrix[0])):
            for k in range(len(second_matrix)):
                result[i][j] += first_matrix[i][k] * second_matrix[k][j]
    return result

def write_output_file(filename, data_type, result_matrix):
    with open(filename, 'w') as file:
        file.write(f"{data_type}\n")
        file.write(f"{len(result_matrix)} {len(result_matrix[0])}\n")
        for row in result_matrix:
            file.write(' '.join(map(str, row)) + '\n')

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py input_filename output_filename")
        return

    filename_in = sys.argv[1]
    filename_out = sys.argv[2]

    opcode, data_type, first_matrix, second_matrix = read_input_file(filename_in)

    start_time = time.time()
    if opcode == 1:
        result = multiply_matrices(first_matrix, second_matrix)
    else:
        print("Opcoode not supported in Python.")
        return

    execution_time = time.time() - start_time

    write_output_file(filename_out, data_type, result)

    print("Output written to:", filename_out)
    print("Execution time:", execution_time, "seconds")

if __name__ == "__main__":
    main()




