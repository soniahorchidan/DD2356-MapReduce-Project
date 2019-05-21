#include <ctype.h>
#include <omp.h>
#include "framework.h"

typedef struct {
	char* key;
	int value;
}KeyValue;

typedef struct {

	// input, output files
	MPI_File input_file, output_file;

	// Datatype for reading/writing
	MPI_Datatype subarray;	

	// local data
	KeyValue *data;

	// communicators
	MPI_Comm grid_comm;
	MPI_Comm row_comm;
	MPI_Comm col_comm;

	// ranks ans sizes
	int world_rank, world_size;
	int grid_rank;
	int grid_dim[2]; // split into n x n processes
	int grid_coords[2];
	int row_rank, row_size;
	int col_rank, col_size;

	// input size
	long input_len;
	//local size
	long local_data_len;
	// number of extra chars to read on both sides
	int offset;


}LocalConfig;

LocalConfig lc;

void read_file(char *input){

	// get world details
	MPI_Comm_rank(MPI_COMM_WORLD, &lc.world_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &lc.world_size);

	//  make sure  world dimensions are square
	if ((int)sqrt(lc.world_size) !=  sqrt(lc.world_size)){
		printf("Processes size is not square\n");
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	lc.grid_dim[0] = lc.grid_dim[1] = sqrt(lc.world_size);

	// read input file size
	if(lc.world_rank == 0) {
	MPI_File_open(MPI_COMM_SELF, input, MPI_MODE_RDONLY, MPI_INFO_NULL, &lc.input_file);
	MPI_Offset size;
	MPI_File_get_size(lc.input_file, &size);
	MPI_File_close(&lc.input_file);
	lc.input_len = size;
	printf("Input file contains %ld chars.\n", lc.input_len);
	}

	// broadcast input size
	MPI_Bcast(&lc.input_len, 1, MPI_LONG, 0, MPI_COMM_WORLD);

	// create communicators for NxN processes
	int period[2] = {1, 1};
	MPI_Cart_create(MPI_COMM_WORLD, 2, lc.grid_dim, period , 1, &lc.grid_comm);

	MPI_Cart_coords(lc.grid_comm, lc.world_rank, 2, lc.grid_coords);
	MPI_Cart_rank(lc.grid_comm, lc.grid_coords, &lc.grid_rank);


	// Sub div cart communicator to N row communicator
	int selected_dim[2] = {0, 1};
	MPI_Cart_sub(lc.grid_comm, selected_dim, &lc.row_comm);
	MPI_Comm_rank(lc.row_comm, &lc.row_rank);
	MPI_Comm_size(lc.row_comm, &lc.row_size);


	// Sub div cart communicator to N col communicator
	selected_dim[0] = 1;
	selected_dim[1] = 0;
	MPI_Cart_sub(lc.grid_comm, selected_dim, &lc.col_comm);
	MPI_Comm_rank(lc.col_comm, &lc.col_rank);
	MPI_Comm_size(lc.col_comm, &lc.col_size);

	// calculate local data size
	lc.offset = 2;
	int chunk_size = lc.input_len / lc.world_size;
	int read_len = chunk_size;
	read_len += 2*lc.offset; // add size of sides
	
	// first process cannot read offset at the beginning of file
	if (lc.world_rank == 0) read_len -= lc.offset;
	
	// last process cannot read offset at the end of file
	// it should also read any remaining bytes
	if (lc.world_rank == lc.world_size -1) {
		read_len -= lc.offset;
		read_len += lc.input_len % lc.world_size;
	}
	
	// create subarray datatype
	int start_from[1] = {lc.world_rank * chunk_size};
	int array_size[1] = {lc.input_len};
	int subarray_size[1] = {read_len};
	if (lc.world_rank != 0) start_from[0]-= lc.offset;
	
	MPI_Type_create_subarray(1, array_size, subarray_size, start_from, MPI_ORDER_C, MPI_CHAR, &lc.subarray);
	MPI_Type_commit(&lc.subarray);

	// alloc space for reading
	char *p = (char*) malloc((read_len+1) * sizeof(char));

	// read file
	MPI_File_open(MPI_COMM_WORLD, input, MPI_MODE_RDONLY, MPI_INFO_NULL, &lc.input_file);
	MPI_File_set_view(lc.input_file, 0, MPI_CHAR, lc.subarray, "native", MPI_INFO_NULL);
	MPI_File_read_all(lc.input_file, p, read_len, MPI_CHAR, MPI_STATUS_IGNORE);
	MPI_File_close(&lc.input_file);

	// terminate string
	p[read_len] = '\0';
	printf("Rank %d read: %s\n",lc.world_rank ,p);

	lc.data = (KeyValue*) malloc(sizeof(KeyValue));
	lc.data[0].key = p;
	lc.data[0].value = 0;
	lc.local_data_len = 1;

}


int isseparator(char c) {
	char *separators = ".,;:\"'()[]{}<>/?!\\\n ";
	return strchr(separators, c) != NULL ? 1 : 0;
}


void find_next_word(char* input_buffer, int *start_index, int *end_index) {
	int current_index = (*start_index);
	int input_length =  strlen(lc.data[0].key);
	
	while(current_index < input_length) {

		while(isseparator(input_buffer[current_index]) && current_index < input_length)
			current_index ++;

		*start_index = current_index;		//start of first word; hopefully

		if(isdigit(input_buffer[current_index])) {	// should be followed only by digits
			while(isdigit(input_buffer[current_index]) && current_index < input_length)
				current_index ++;
			if(isseparator(input_buffer[current_index])) {	// found word
				*end_index = current_index - 1;
				current_index = input_length;
			} 
			else { // not a valid word; skip until separator
				while(!isseparator(input_buffer[current_index]) && current_index < input_length)
					current_index ++;
			}
		} else 	if(isalpha(input_buffer[current_index])) {	// should be followed only by digits
			while(isalpha(input_buffer[current_index]) && current_index < input_length)
				current_index ++;
			if(isseparator(input_buffer[current_index])) {	// found word
				*end_index = current_index - 1;
				current_index = input_length;
			} 
			else { // not a valid word; skip until separator
				while(!isseparator(input_buffer[current_index]) && current_index < input_length)
					current_index ++;
			}
		}
	}
	
}


void flat_map(){

	int start_index = 0;
	int end_index = 0;
	char **words = (char **)calloc(100, sizeof(char *));
	int i;
	for (i = 0; i < 100; i ++)
   	     words[i] = (char *)calloc(100, sizeof(char));

	int index = 0;

	long int my_length = strlen(lc.data[0].key);

	while(start_index < my_length - 1) {
		find_next_word(lc.data[0].key, &start_index, &end_index);
		int word_size = end_index - start_index + 1;
		if(word_size <= 0) break;
		words[index] = (char *)realloc(words[index], (word_size + 2) * sizeof(char));
		strncpy(words[index], lc.data[0].key + start_index, word_size);
		words[index][word_size] = '\0';
		index ++;
		if (index % 99 == 0)
			words = (char **)realloc(words, sizeof(char **) * index * 3/2);
		start_index = end_index + 1;
	}


	lc.data = (KeyValue*) realloc(lc.data, index * sizeof(KeyValue));

	#pragma omp parallel for private(i)
	for(i = 0; i < index; i ++) {
		KeyValue new_kv_pair;
		new_kv_pair.key = (char *)malloc(strlen(words[i]) * sizeof(char));
		strcpy(new_kv_pair.key, words[i]);
		new_kv_pair.value = 1;
		lc.data[i] = new_kv_pair;
	}

	lc.local_data_len = index - 1;

	for(i = 0; i < index; i ++) 
		printf("%d: %s %d\n", lc.world_rank, lc.data[i].key, lc.data[i].value);
	
}

void reduce(){
	
}

void write_file(){

}


