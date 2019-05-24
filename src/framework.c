#include <ctype.h>
#include <omp.h>
#include "framework.h"

typedef struct {
    char *key;
    int value;
} KeyValue;

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


} LocalConfig;

LocalConfig lc;

void read_file(char *input) {

    // get world details
    MPI_Comm_rank(MPI_COMM_WORLD, &lc.world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &lc.world_size);

    //  make sure  world dimensions are square
    if ((int) sqrt(lc.world_size) != sqrt(lc.world_size)) {
        printf("Processes size is not square\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    lc.grid_dim[0] = lc.grid_dim[1] = sqrt(lc.world_size);

    // read input file size
    if (lc.world_rank == 0) {
        MPI_File_open(MPI_COMM_SELF, input, MPI_MODE_RDONLY, MPI_INFO_NULL, &lc.input_file);
        MPI_Offset size;
        MPI_File_get_size(lc.input_file, &size);
        MPI_File_close(&lc.input_file);
        lc.input_len = size - 1; // exclude EOF
        printf("Input file contains %ld chars.\n", lc.input_len);
    }

    // broadcast input size
    MPI_Bcast(&lc.input_len, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    // create communicators for NxN processes
    int period[2] = {1, 1};
    MPI_Cart_create(MPI_COMM_WORLD, 2, lc.grid_dim, period, 1, &lc.grid_comm);

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
    read_len += 2 * lc.offset; // add size of sides

    // first process cannot read offset at the beginning of file
    if (lc.world_rank == 0) read_len -= lc.offset;

    // last process cannot read offset at the end of file
    // it should also read any remaining bytes
    if (lc.world_rank == lc.world_size - 1) {
        read_len -= lc.offset;
        read_len += lc.input_len % lc.world_size;
    }
    // create subarray datatype
    int start_from[1] = {lc.world_rank * chunk_size};
    int array_size[1] = {lc.input_len};
    int subarray_size[1] = {read_len};
    if (lc.world_rank != 0) start_from[0] -= lc.offset;

    MPI_Type_create_subarray(1, array_size, subarray_size, start_from, MPI_ORDER_C, MPI_CHAR, &lc.subarray);
    MPI_Type_commit(&lc.subarray);

    // alloc space for reading
    char *p = (char *) malloc((read_len) * sizeof(char));

    // read file
    MPI_File_open(MPI_COMM_WORLD, input, MPI_MODE_RDONLY, MPI_INFO_NULL, &lc.input_file);
    MPI_File_set_view(lc.input_file, 0, MPI_CHAR, lc.subarray, "native", MPI_INFO_NULL);
    MPI_File_read_all(lc.input_file, p, read_len, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&lc.input_file);

    //printf("Rank %d read: |%s|\n",lc.world_rank ,p);

    lc.data = (KeyValue *) malloc(sizeof(KeyValue));
    lc.data[0].key = p;
    lc.data[0].value = 1;
    lc.local_data_len = 1;

}


int isseparator(char c) {
    //char *separators = ".,;:\"'()[]{}<>/?!\"\\\n ";
    // return strchr(separators, c) != NULL ? 1 : 0;
    return !isdigit(c) && !isalpha(c);
}


void find_next_word(char *input_buffer, int *start_index, int *end_index) {
    int current_index = (*start_index);
    int input_length = strlen(lc.data[0].key);

    while (current_index < input_length) {

        while (isseparator(input_buffer[current_index]) && current_index < input_length)
            current_index++;

        *start_index = current_index;        //start of first word; hopefully

        if (isdigit(input_buffer[current_index])) {    // should be followed only by digits
            while (isdigit(input_buffer[current_index]) && current_index < input_length)
                current_index++;
            if (isseparator(input_buffer[current_index])) {    // found word
                *end_index = current_index - 1;
                current_index = input_length;
            } else { // not a valid word; skip until separator
                while (!isseparator(input_buffer[current_index]) && current_index < input_length)
                    current_index++;
            }
        } else if (isalpha(input_buffer[current_index])) {    // should be followed only by digits
            while (isalpha(input_buffer[current_index]) && current_index < input_length)
                current_index++;
            if (isseparator(input_buffer[current_index])) {    // found word
                *end_index = current_index - 1;
                current_index = input_length;
            } else { // not a valid word; skip until separator
                while (!isseparator(input_buffer[current_index]) && current_index < input_length)
                    current_index++;
            }
        }
    }

}


void flat_map() {

    int start_index = 0;
    int end_index = 0;
    char **words = (char **) malloc(100 * sizeof(char *));
    int i;
    int index = 0;

    long int my_length = strlen(lc.data[0].key);
    while (start_index < my_length - 1) {
        find_next_word(lc.data[0].key, &start_index, &end_index);
        int word_size = end_index - start_index + 1;
        if (word_size <= 0) break;
        words[index] = (char *) malloc((word_size + 1) * sizeof(char));
        strncpy(words[index], lc.data[0].key + start_index, word_size);
        words[index][word_size] = '\0';
        index++;
        if (index % 90 == 0)
            words = (char **) realloc(words, 2 * index * sizeof(*words));
        start_index = end_index + 1;
    }


    lc.data = (KeyValue *) realloc(lc.data, index * sizeof(KeyValue));

#pragma omp parallel for private(i)
    for (i = 0; i < index; i++) {
        KeyValue new_kv_pair;
        new_kv_pair.key = (char *) malloc(strlen(words[i]) * sizeof(char));
        strcpy(new_kv_pair.key, words[i]);
        new_kv_pair.value = 1;
        lc.data[i] = new_kv_pair;
    }

    lc.local_data_len = index;

    // for(i = 0; i < index; i ++)
    // 	printf("%d: %s %d\n", lc.world_rank, lc.data[i].key, lc.data[i].value);

    free(words);
}

unsigned long hash(char *str) {
    unsigned long hash = 5381;
    int c;

    while ((c = *str++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}

void reduce_local() {

    if (lc.local_data_len == 1) return; // nothing to merge

    int merged = 0; //count how many elements were reduced
    int i, j;

    for (i = lc.local_data_len - 1; i > 0; i--) {

        int s1 = strlen(lc.data[i].key);

        // compare with others
        for (j = 0; j < i; j++) {

            int s2 = strlen(lc.data[j].key);

            // if keys are equal
            if (s1 == s2 && memcmp(lc.data[i].key, lc.data[j].key, s1) == 0) {
                free(lc.data[i].key);
                lc.data[i].key = NULL;
                lc.data[j].value++;
                merged++;
                break;
            }
        }
    }

    // copy to new bucket
    j = 0;
    KeyValue *new_data = (KeyValue *) malloc((lc.local_data_len - merged) * sizeof(KeyValue));
    for (i = 0; i < lc.local_data_len; i++) {
        if (lc.data[i].key != NULL) {
            if (lc.world_rank == 0) printf("NOT NULL\n");
            new_data[j].key = lc.data[i].key;
            new_data[j].value = lc.data[i].value;
            j++;
        }
    }

    free(lc.data);
    lc.data = new_data;
    lc.local_data_len -= merged;
}

void reduce() {

    // local reduce
    local_reduce();

    KeyValue **buckets = (KeyValue **) malloc(lc.world_size * sizeof(KeyValue *));
    int i;
    for (i = 0; i < lc.world_size; i++)
        buckets[i] = (KeyValue *) malloc(sizeof(KeyValue));

    int *sizes = (int *) calloc(lc.world_size, sizeof(int));
    for (i = 0; i < lc.local_data_len; i++) {
        unsigned long word_hash = hash(lc.data[i].key);
        int receiver = word_hash % lc.world_size;
        buckets[receiver][sizes[receiver]++] = lc.data[i];
        buckets[receiver] = (KeyValue *) realloc(buckets[receiver], (sizes[receiver] + 1) * sizeof(KeyValue));
    }

    // printf("proc %d: ", lc.world_rank);
    // for(i = 0; i < lc.world_size; i ++)
    //     printf("%d ", sizes[i]);
    // printf("\n\n");


    int *rec_sizes = (int *) calloc(lc.world_size, sizeof(int));

    MPI_Alltoall(sizes, 1, MPI_INT, rec_sizes, 1, MPI_INT, MPI_COMM_WORLD);

    // printf("proc %d: ", lc.world_rank);
    // for(i = 0; i < lc.world_size; i ++)
    //     printf("%d ", rec_sizes[i]);
    // printf("\n\n");

    free(sizes);
    free(buckets);

}


void write_file() {
    //calculate local size
    int local_size = 0;
    int i;
    for (i = 0; i < lc.local_data_len; i++) {
        KeyValue kv = lc.data[i];
        // size of string
        local_size += strlen(kv.key);
        // size of count
        if (kv.value == 0) {
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
        local_size += floor(log10(abs(kv.value))) + 1;
        // space and newline
        local_size += 2;
    }
    // printf("Rank %d will print %d chars\n", lc.world_rank, local_size);

    // create local result
    char *result = (char *) malloc(local_size * sizeof(char));
    int j = 0;
    for (i = 0; i < lc.local_data_len; i++) {
        // get key value pair
        KeyValue kv = lc.data[i];
        j += sprintf(&result[j], "%s %d\n", kv.key, kv.value);

    }
    // printf("Rank %d result: |%s|\n", lc.world_rank, result);


    // local sizes are distributed across the network
    int proc_size[lc.world_size];
    MPI_Allgather(&local_size, 1, MPI_INT, &proc_size, 1, MPI_INT, MPI_COMM_WORLD);

    // calculate offset
    int proc_offset = 0;
    int total_size = 0;
    i = 0;
    while (i < lc.world_size) {
        if (i < lc.world_rank) proc_offset += proc_size[i];
        total_size += proc_size[i];
        i++;
    }
    //printf("Rank %d offset %d chars\n", lc.world_rank, proc_offset);

    // write results to file
    MPI_Datatype result_datatype;
    int start_from[1] = {proc_offset};
    int local_array[1] = {local_size};
    int result_array[1] = {total_size};
    MPI_Type_create_subarray(1, result_array, local_array, start_from, MPI_ORDER_C, MPI_CHAR, &result_datatype);
    MPI_Type_commit(&result_datatype);
    MPI_File_open(MPI_COMM_WORLD, "result", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &lc.output_file);
    MPI_File_set_view(lc.output_file, 0, MPI_CHAR, result_datatype, "native", MPI_INFO_NULL);
    MPI_File_write_all(lc.output_file, result, local_size, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close(&lc.output_file);

    // clean up
    free(result);
    free(lc.data);

}


