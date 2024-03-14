#include <stdio.h>
#include <stdlib.h>
int main(int argc,char *argv[])
{
	int opcode=0;
	FILE *inputFile = fopen(argv[1], "r");
	if (inputFile == NULL) {
    	perror("Error opening input file");
    	return 1;
	}
	fscanf(inputFile, "%d", &opcode);
	if(opcode==1)
    	{
   		 printf("Running Python code\n");
   		 char buffer[100];
   		 char* command ="python3 matrix_mul.py";
   		 int j=snprintf(buffer, 99, "%s %s %s", command,argv[1],argv[2] );
   		 system(buffer);
    	}
   if(opcode==2)
    	{
   		 printf("Running C code\n");
   		 char buffer[100];
   		 char* command ="./cmm ";
   		 int j=snprintf(buffer, 99, "%s %s %s", command,argv[1],argv[2] );
   		 system(buffer);
    	}
   if(opcode==3)
    	{
   		 printf("Running Java code\n");
   		 char buffer[100];
   		 char* command ="java MatrixMultiplication ";
   		 int j=snprintf(buffer, 99, "%s %s %s", command,argv[1],argv[2] );
   		 system(buffer);
    	}
	if(opcode==4)
    	{
   		 printf("Running C code with vector instructions\n");
   		 char buffer[100];
   		 char* command ="./mv ";
   		 int j=snprintf(buffer, 99, "%s %s %s", command,argv[1],argv[2] );
   		 system(buffer);
    	}
	if(opcode==5)
    	{
   		 printf("Running C code with Pthreads\n");
   		 char buffer[100];
   		 char* command ="./mmt ";
   		 int j=snprintf(buffer, 99, "%s %s %s", command,argv[1],argv[2] );
   		 system(buffer);
    	}
	if(opcode==6)
    	{
   		 printf("Running C code with vector instructions+Pthreads\n");
   		 char buffer[100];
   		 char* command ="./mvt ";
   		 int j=snprintf(buffer, 99, "%s %s %s", command,argv[1],argv[2] );
   		 system(buffer);
    	}
   	 
	else
	{
	printf("Opcode not supported\n");
	}
   printf("Thanks For using\n");


}

