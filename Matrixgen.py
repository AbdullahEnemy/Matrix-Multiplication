import random
import numpy as np




def generate_matrix(rows, cols, data_type):
	if data_type == int:
		return np.random.randint(-9223372036854775808, 9223372036854775808, size=(rows, cols), dtype=np.int64)

	else:
		return np.random.uniform(-2**63,2**63 - 1, size=(rows, cols)).astype(np.float64)
	
def generate_matrix_int32(rows, cols):
	return np.random.randint(-2**31, 2**31 - 1, size=(rows, cols), dtype=np.int32)


def generate_matrix_float32(rows, cols):
	return np.random.uniform(-2**31,2**31 - 1, size=(rows, cols)).astype(np.float32)


def write_matrix_to_file(op_code, matrix1, matrix2, data_type, file):
	with open(file, 'w') as f:
		f.write(str(op_code) + '\n')
		f.write(str(data_type) + '\n')
		f.write(f"{len(matrix1)} {len(matrix1[0])}\n")
		for row in matrix1:
			f.write(" ".join(map(str, row)) + '\n')
		
		f.write(f"{len(matrix2)} {len(matrix2[0])}\n")
		for row in matrix2:
			f.write(" ".join(map(str, row)) + '\n')
    	
        	





print('Enter Required Data type. i.e., 1.int, 2.float, 3.longint, 4.double')
typedata = int(input())
if typedata == 1:
	data_type = int
elif typedata == 2:
	data_type = float
elif typedata == 3:
	data_type = int
elif typedata == 4:
	data_type = float




print('Enter Required Opcode. i.e., 1.Python, 2.C, 3.Java, 4.C vector instructions, 5.pthreads based, 6.with pthreads based multithreading and vector instructions together')
opcode = int(input())
print('Enter dimensions for matrix1:')
rows1, cols1 = map(int, input().split())
print('Enter dimensions for matrix2:')
rows2, cols2 = map(int, input().split())
# Check data type and generate matrices accordingly
if typedata in [3, 4]:
	matrix1 = generate_matrix(rows1, cols1, data_type)
	matrix2 = generate_matrix(rows2, cols2, data_type)
elif typedata == 1:
	matrix1 = generate_matrix_int32(rows1, cols1)
	matrix2 = generate_matrix_int32(rows2, cols2)
elif typedata == 2:
	matrix1 = generate_matrix_float32(rows1, cols1)
	matrix2 = generate_matrix_float32(rows2, cols2)


print('Enter output FileName:')
filename = input()
write_matrix_to_file(opcode, matrix1, matrix2, typedata, filename)







