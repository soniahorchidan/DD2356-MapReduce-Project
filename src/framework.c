#include "framework.h"
#include "structs.h"

struct Data {
	/* MPI Files */
	MPI_File A_file, B_file, C_file;
	char *outfile;

	/* MPI Datatypes for matrix blocks */
	MPI_Datatype block;

	/* Matrix data */
	double *A, *A_tmp, *B, *B_tmp, *C;

	/* Cart communicators */
	MPI_Comm grid_comm;
	MPI_Comm row_comm;
	MPI_Comm col_comm;

	/* Cart communicator dim and ranks */
	int dim[2], coords[2];
	int world_rank, world_size, grid_rank;
	int row_rank, row_size, col_rank, col_size;

	/* Full matrix dim */
	int A_dims[2];
	int B_dims[2];
	int C_dims[2];
	int matrix_size;

	/* Process local matrix dim */
	int local_dims[2];
	int local_size;
};

struct Data data;

int read_file(char *input){
	
}

int flat_map(int data, Chuck* (*mapper)(Chuck)){
	
}

int reduce(int data, Chuck (*reducer)(Chuck, Chunk)){
	
}

int write_file(int data, char *output){
	
}

void free(int data){
	
}


