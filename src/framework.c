#include <ctype.h>

#include <omp.h>

#include "framework.h"

typedef struct {
    char * key;
    int value;
}
KeyValue;

typedef struct {

    // input, output files
    MPI_File input_file, output_file;

    // Datatype for reading/writing
    MPI_Datatype subarray;

    // local data
    KeyValue * data;
    KeyValue input_data;

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
    // number of extra chars to read in the end of array
    int offset;
    // max word size
    int max_word_size;

}
LocalConfig;

LocalConfig lc;

void read_file(char * input) {

     lc.max_word_size = 16;

    // get world details
    MPI_Comm_rank(MPI_COMM_WORLD, & lc.world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, & lc.world_size);

    //  make sure  world dimensions are square
    if ((int) sqrt(lc.world_size) != sqrt(lc.world_size)) {
        printf("Processes size is not square\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    lc.grid_dim[0] = lc.grid_dim[1] = sqrt(lc.world_size);

    // read input file size
    if (lc.world_rank == 0) {
        MPI_File_open(MPI_COMM_SELF, input, MPI_MODE_RDONLY, MPI_INFO_NULL, & lc.input_file);
        MPI_Offset size;
        MPI_File_get_size(lc.input_file, & size);
        MPI_File_close( & lc.input_file);
        lc.input_len = size - 1; // exclude EOF
        // printf("Input file contains %ld chars.\n", lc.input_len);
    }

    // broadcast input size
    MPI_Bcast( & lc.input_len, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    // create communicators for NxN processes
    int period[2] = {
        1,
        1
    };

    MPI_Cart_create(MPI_COMM_WORLD, 2, lc.grid_dim, period, 1, & lc.grid_comm);

    MPI_Cart_coords(lc.grid_comm, lc.world_rank, 2, lc.grid_coords);
    MPI_Cart_rank(lc.grid_comm, lc.grid_coords, & lc.grid_rank);

    // Sub div cart communicator to N row communicator
    int selected_dim[2] = {
        0,
        1
    };
    MPI_Cart_sub(lc.grid_comm, selected_dim, & lc.row_comm);
    MPI_Comm_rank(lc.row_comm, & lc.row_rank);
    MPI_Comm_size(lc.row_comm, & lc.row_size);

    // Sub div cart communicator to N col communicator
    selected_dim[0] = 1;
    selected_dim[1] = 0;
    MPI_Cart_sub(lc.grid_comm, selected_dim, & lc.col_comm);
    MPI_Comm_rank(lc.col_comm, & lc.col_rank);
    MPI_Comm_size(lc.col_comm, & lc.col_size);

    // calculate local data size
    lc.offset = 15;
    int chunk_size = 64; //lc.input_len / lc.world_size;
    int read_len = chunk_size;
    read_len += lc.offset + 1; // +1 for the beginning

    // first process cannot read offset at the beginning of file
    if (lc.world_rank == 0) read_len -= 1;

    // last process cannot read offset at the end of file
    // it should also read any remaining bytes
    if (lc.world_rank == lc.world_size - 1) {
        read_len -= lc.offset;
        read_len += lc.input_len % lc.world_size;
    }

    // alloc space for reading
    char *p, *p_dummy;

    if(lc.world_rank == 0){
        p_dummy = (char * ) malloc((read_len+2) * sizeof(char));
        p = &p_dummy[1]; // we save the first slot for a separator
	p_dummy[0] = ' ';
        p_dummy[read_len+1] = '\0';
    }else{
        p = (char * ) malloc((read_len+1) * sizeof(char));
        p[read_len] = '\0';
    }

    long num_reads = lc.input_len / (chunk_size * lc.world_size);
    long remaining = lc.input_len  - num_reads * chunk_size * lc.world_size;
    long extra_reads = remaining / chunk_size;
    if(lc.world_rank <= extra_reads) {
        num_reads ++;
    }

    MPI_Aint length = chunk_size * sizeof(char);
    MPI_Offset disp = lc.world_rank * length;
    MPI_Aint extent = lc.world_size * length;
    MPI_Datatype contig, filetype;

    MPI_Type_contiguous(read_len, MPI_CHAR, &contig);
    MPI_Type_create_resized(contig, 0, extent, &filetype);
    MPI_Type_commit(&filetype);

    // read file
    MPI_File_open(MPI_COMM_WORLD, input, MPI_MODE_RDONLY, MPI_INFO_NULL, & lc.input_file);
    MPI_File_set_view(lc.input_file, disp, MPI_CHAR, filetype, "native", MPI_INFO_NULL);

    lc.data = (KeyValue * ) malloc(sizeof(KeyValue));
    lc.local_data_len = 0;

    int i;
    
    for(i = 0; i < num_reads; i ++) {
          if(i == num_reads - 1) read_len = remaining;
          MPI_File_read(lc.input_file, p, read_len, MPI_CHAR, MPI_STATUS_IGNORE);
          if(lc.world_rank == 0) p = p_dummy;
          // printf("Rank %d read: |%s|\n",lc.world_rank ,p);
              
          lc.input_data.key = p;
          lc.input_data.value = 1;
          flat_map();
    }

    MPI_File_close( & lc.input_file);

}

int isSep(char c) {
    //char *separators = ".,;:\"'()[]{}<>/?!\"\\\n ";
    // return strchr(separators, c) != NULL ? 1 : 0;
    return !isdigit(c) && !isalpha(c);
}

void find_next_word(char * input_buffer, int * read_head, int * word_start_index, int * word_size) {
    int total_chunk_len = strlen(input_buffer);
    int chunk_len = total_chunk_len - lc.offset;
    if (lc.world_rank == lc.world_size -1) chunk_len+= lc.offset;

    while (isSep(input_buffer[*read_head]) && *read_head < chunk_len) (*read_head)++;

    *word_size = 0;
    // if we are in next chunk's region
    if(*read_head == chunk_len) return;

    *word_start_index = *read_head; //start of first word; hopefully
    do{
        (*word_size)++;
        (*read_head)++;
    }while( (*word_size < lc.max_word_size) && !isSep(input_buffer[*read_head]) && (*read_head < total_chunk_len));
    while(!isSep(input_buffer[*read_head]) && *read_head < chunk_len) (*read_head)++;
}



void flat_map() {

    int buffer_size = 100;
    int word_start_index;
    int word_size;
    char ** words = (char ** ) malloc(buffer_size * sizeof(char * ));
    int i;
    int word_counter = 0;
    int read_head = 0;

    long int chunk_len = strlen(lc.input_data.key) - lc.offset; // remove offset size
    if(lc.world_rank == lc.world_size -1) chunk_len+= lc.offset; // last chuck has no offset

    // skip broken word - if any
    while(!isSep(lc.input_data.key[read_head]) && read_head < chunk_len) read_head++;

    while (read_head < chunk_len) {
        find_next_word(lc.input_data.key, &read_head, & word_start_index, & word_size);
        if (word_size == 0) break;
        words[word_counter] = (char * ) malloc((word_size + 1) * sizeof(char));
        strncpy(words[word_counter], lc.input_data.key + word_start_index, word_size);
        words[word_counter][word_size] = '\0';
        word_counter++;
        if (word_counter % buffer_size == 0) // full 
            words = (char ** ) realloc(words, 2 * word_counter * sizeof( * words));
    }

    lc.data = (KeyValue * ) realloc(lc.data, (word_counter + lc.local_data_len + 1) * sizeof(KeyValue));

    #pragma omp parallel for private(i)
    for (i = lc.local_data_len; i < lc.local_data_len + word_counter; i++) {
        KeyValue new_kv_pair;
        new_kv_pair.key = (char * ) malloc(strlen(words[i - lc.local_data_len]) * sizeof(char));
        strcpy(new_kv_pair.key, words[i - lc.local_data_len]);
        new_kv_pair.value = 1;
        lc.data[i] = new_kv_pair;
    }

    lc.local_data_len += word_counter;
    free(words);
}

unsigned long hash(char * str) {
    unsigned long hash = 5381;
    int c;

    while ((c = * str++))
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
    KeyValue * new_data = (KeyValue * ) malloc((lc.local_data_len - merged) * sizeof(KeyValue));
    for (i = 0; i < lc.local_data_len; i++) {
        if (lc.data[i].key != NULL) {
            new_data[j].key = lc.data[i].key;
            new_data[j].value = lc.data[i].value;
            j++;
        }
    }

    free(lc.data);
    lc.data = new_data;
    lc.local_data_len -= merged;
}

void string_to_byte_array(char * input, char * output) {
    int loop = 0;
    int i = 0;
    while (input[loop] != '\0') {
        output[i++] = input[loop++];
    }
}

void int_to_byte_array(int input, char * output) {
    output[3] = (input >> 24) & 0xFF;
    output[2] = (input >> 16) & 0xFF;
    output[1] = (input >> 8) & 0xFF;
    output[0] = input & 0xFF;
}

void convert_buckets_into_bytes(KeyValue * bucket, int bucket_size, char ** bytes, int * bytes_size) {
    int i;
    * bytes_size = 0;

    for (i = 0; i < bucket_size; i++) {
        char * key = (char * ) malloc(strlen(bucket[i].key) * sizeof(char));
        string_to_byte_array(bucket[i].key, key);

        char * value = (char * ) malloc(sizeof(bucket[i].value) * sizeof(char));
        int_to_byte_array(bucket[i].value, value);

        int current_size = strlen(bucket[i].key) + sizeof(bucket[i].value) + 1;

        char b[current_size * sizeof(char) + 1];
        b[0] = '\0';
        strcat(b, key);

        b[strlen(bucket[i].key)] = '\0';
        memcpy(b + strlen(bucket[i].key) + 1, value, sizeof(int));

        free(key);
        free(value);

        * bytes = (char * ) realloc( * bytes, ( * bytes_size + current_size) * sizeof(char));
        memcpy( * bytes + * bytes_size, b, current_size);
        * bytes_size += current_size;

    }
}

void merge(char * recv, int recv_size) {

    char buffer[16];
    int j, i, k;
    int remaining = 0;

    i = 0;
    while (i < recv_size) {

        // read key
        j = 0;
        while (recv[i] != '\0') {
            buffer[j] = recv[i];
            i++;
            j++;
        }
        i++; // for '\0'
        buffer[j] = '\0';

        // read value
        int value = recv[i] | recv[i + 1] << 8 | recv[i + 2] << 16 | recv[i + 3] << 24;

        //if (lc.world_rank == 0)printf("entry: |%s|, %d\n", buffer, value);

        // try to merge with local bucket
        int merged = 0;
        for (k = 0; k < lc.local_data_len; k++) {
            if (strlen(lc.data[k].key) == j && memcmp(lc.data[k].key, buffer, j) == 0) {
                // merge
                lc.data[k].value += value;
                // set value to -1
                recv[i] = -1;
                recv[i + 1] = recv[i + 2] = recv[i + 3] = 0;

                merged = 1;
                break;
            }
        }

        if (merged == 0) {
            remaining++;
        }
        i += 4;
    }

    if (remaining == 0) return; // nothing to merge
    int merge_size = lc.local_data_len + remaining;
    KeyValue * merged = (KeyValue * ) malloc(merge_size * sizeof(KeyValue));

    // copy bucket
    for (i = 0; i < lc.local_data_len; i++) {
        merged[i].key = lc.data[i].key;
        merged[i].value = lc.data[i].value;
    }

    // append rest of words to bucket
    i = 0;
    int l = lc.local_data_len;
    while (i < recv_size) {

        // read key
        j = 0;
        while (recv[i] != '\0') {
            buffer[j] = recv[i];
            i++;
            j++;
        }
        i++; // for '\0'
        buffer[j] = '\0';
        j++;

        // read value
        int value = recv[i] | recv[i + 1] << 8 | recv[i + 2] << 16 | recv[i + 3] << 24;
        if (value != -1) {
            char * array = (char * ) malloc(j);
            memcpy(array, buffer, j);
            merged[l].key = array;
            merged[l].value = value;
            l++;
        }
        i += 4;
    }

    free(lc.data);
    lc.data = merged;
    lc.local_data_len = merge_size;

}

void reduce() {

    // local reduce
    reduce_local();

    KeyValue ** buckets = (KeyValue ** ) malloc(lc.world_size * sizeof(KeyValue * ));
    int i;
    for (i = 0; i < lc.world_size; i++)
        buckets[i] = (KeyValue * ) malloc(sizeof(KeyValue));

    // find buckets
    int * sizes = (int * ) calloc(lc.world_size, sizeof(int));
    for (i = 0; i < lc.local_data_len; i++) {
        unsigned long word_hash = hash(lc.data[i].key);
        int receiver = word_hash % lc.world_size;
        buckets[receiver][sizes[receiver]++] = lc.data[i];
        buckets[receiver] = (KeyValue * ) realloc(buckets[receiver], (sizes[receiver] + 1) * sizeof(KeyValue));
    }

    int * sizes_bytes = (int * ) calloc(lc.world_size, sizeof(int));
    char ** send_bytes = (char ** ) malloc(lc.world_size * sizeof(char * ));

    // convert buckets to byte arrays
    for (i = 0; i < lc.world_size; i++) {
        send_bytes[i] = NULL;
        convert_buckets_into_bytes(buckets[i], sizes[i], & send_bytes[i], & sizes_bytes[i]);
    }


    int * recv_sizes_bytes = (int * ) calloc(lc.world_size, sizeof(int));

    // exchange sizes between processes
    MPI_Alltoall(sizes_bytes, 1, MPI_INT, recv_sizes_bytes, 1, MPI_INT, MPI_COMM_WORLD);

     
    MPI_Request send_requests[lc.world_size];

    for (i = 0; i < lc.world_size; i++) {
        // find the right communicator to send to
        if(lc.world_rank % lc.grid_dim[0] == i % lc.grid_dim[0])    //same col
            MPI_Isend(send_bytes[i], sizes_bytes[i], MPI_BYTE, i / lc.grid_dim[0], 0, lc.col_comm, &send_requests[i]);
        else if(lc.world_rank / lc.grid_dim[0] == i / lc.grid_dim[0])    //same row
            MPI_Isend(send_bytes[i], sizes_bytes[i], MPI_BYTE, i % lc.grid_dim[0], 0, lc.row_comm, &send_requests[i]);
        else 
            MPI_Isend(send_bytes[i], sizes_bytes[i], MPI_BYTE, i, 0, lc.grid_comm, &send_requests[i]);
    }

    char ** recv_bytes = (char ** ) malloc(lc.world_size * sizeof(char * ));
    MPI_Request recv_requests[lc.world_size];

    for (i = 0; i < lc.world_size; i++) {
        recv_bytes[i] = (char * ) malloc(recv_sizes_bytes[i] * sizeof(char));
        // find the right communicator to receive from
        if(lc.world_rank % lc.grid_dim[0] == i % lc.grid_dim[0])    //same col
            MPI_Irecv(recv_bytes[i], recv_sizes_bytes[i], MPI_BYTE, i / lc.grid_dim[0], 0, lc.col_comm, &recv_requests[i]);
        else if(lc.world_rank / lc.grid_dim[0] == i / lc.grid_dim[0])    //same row
            MPI_Irecv(recv_bytes[i], recv_sizes_bytes[i], MPI_BYTE, i % lc.grid_dim[0], 0, lc.row_comm, &recv_requests[i]);
        else 
            MPI_Irecv(recv_bytes[i], recv_sizes_bytes[i], MPI_BYTE, i, 0, lc.grid_comm, &recv_requests[i]);
    }

    // free(lc.data);
    // lc.local_data_len = 0;
    int index_rec = 0;
    for (i = 0; i < lc.world_size; i ++) {
        // wait for any message to arrive, then merge
        MPI_Waitany(lc.world_size, recv_requests, &index_rec, MPI_STATUS_IGNORE);
        // printf("PROCESS %d RECEIVED FROM PROCESS %d\n", lc.world_rank, index_rec);
        if(index_rec != lc.world_rank) merge(recv_bytes[index_rec], recv_sizes_bytes[index_rec]);
    }

    MPI_Waitall(lc.world_size, send_requests, MPI_STATUS_IGNORE);
 
    free(recv_bytes);
    free(recv_sizes_bytes);
    free(send_bytes);
    free(sizes_bytes);
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
    char * result = (char * ) malloc(local_size * sizeof(char));
    int j = 0;
    for (i = 0; i < lc.local_data_len; i++) {
        // get key value pair
        KeyValue kv = lc.data[i];
        j += sprintf( & result[j], "%s %d\n", kv.key, kv.value);

    }
    // printf("Rank %d result: %s\n\n", lc.world_rank, result);

    // local sizes are distributed across the network
    int proc_size[lc.world_size];
    MPI_Allgather( & local_size, 1, MPI_INT, & proc_size, 1, MPI_INT, MPI_COMM_WORLD);

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
    int start_from[1] = {
        proc_offset
    };
    int local_array[1] = {
        local_size
    };
    int result_array[1] = {
        total_size
    };
    MPI_Type_create_subarray(1, result_array, local_array, start_from, MPI_ORDER_C, MPI_CHAR, & result_datatype);
    MPI_Type_commit( & result_datatype);
    MPI_File_open(MPI_COMM_WORLD, "result.txt", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, & lc.output_file);
    MPI_File_set_view(lc.output_file, 0, MPI_CHAR, result_datatype, "native", MPI_INFO_NULL);
    MPI_File_write_all(lc.output_file, result, local_size, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close( & lc.output_file);

    // clean up
    free(result);
    free(lc.data);
}

