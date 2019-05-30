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

    // ranks ans sizes
    int world_rank, world_size;

    // input size
    long input_len;
    //local size
    long local_data_len;
    // number of extra chars to read in the end of array
    int offset;
    // max word size
    int max_word_size;
    // chunk size per iteration
    int iter_size;

}
LocalConfig;

LocalConfig lc = {.iter_size = 100, .max_word_size = 16};

void read_file(char * input) {

    // get world details
    MPI_Comm_rank(MPI_COMM_WORLD, & lc.world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, & lc.world_size);

    // read input file size
    if (lc.world_rank == 0) {
        MPI_File_open(MPI_COMM_SELF, input, MPI_MODE_RDONLY, MPI_INFO_NULL, & lc.input_file);
        MPI_Offset size;
        MPI_File_get_size(lc.input_file, & size);
        MPI_File_close( & lc.input_file);
        lc.input_len = size - 1; // exclude EOF
        printf("Input file contains %ld chars.\n", lc.input_len);
    }

    // broadcast input size
    MPI_Bcast( & lc.input_len, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    // calculate local data size
    lc.offset = lc.max_word_size -1;
    int chunk_size = lc.input_len / lc.world_size;
    int read_len = chunk_size;

    // last process reads trailing chars - size: [0, lc.world_size)
    if(lc.world_rank == lc.world_size -1) read_len += lc.input_len % lc.world_size;

    //all processes except the first one read one extra byte in the beginning
    if(lc.world_rank != 0) read_len += 1;

    // all processes except the last one read max_word_size-1 chars in the end

    if (lc.world_rank != lc.world_size - 1) {
        read_len += lc.offset;
    }

    // create subarray datatype
    int array_size[2] = {
        lc.world_size+1, chunk_size
    };
    int start_from[2] = {
        lc.world_rank==0? lc.world_rank : lc.world_rank-1, 0
    };
    int subarray_size[2] = {
        lc.world_rank==0? 2 : 3, chunk_size
    };

    MPI_Type_create_subarray(2, array_size, subarray_size, start_from, MPI_ORDER_C, MPI_CHAR, & lc.subarray);
    MPI_Type_commit( & lc.subarray);

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

    int offset = lc.world_rank == 0? 0: chunk_size -1;
    // read file
    MPI_File_open(MPI_COMM_WORLD, input, MPI_MODE_RDONLY, MPI_INFO_NULL, & lc.input_file);
    MPI_File_set_view(lc.input_file, offset, MPI_CHAR, lc.subarray, "native", MPI_INFO_NULL);
    MPI_File_read_all(lc.input_file, p, read_len, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close( & lc.input_file);

    if(lc.world_rank == 0) p = p_dummy;

    //printf("Rank %d read: |%s|\n",lc.world_rank ,p);

    lc.data = (KeyValue * ) malloc(sizeof(KeyValue));
    lc.data[0].key = p;
    lc.data[0].value = 1;
    lc.local_data_len = 1;

}

int isSep(char c) {
    return !isdigit(c) && !isalpha(c);
}

void find_next_word(char * input_buffer, int input_size, int * read_head, int * word_start_index, int * word_size) {
    int chunk_len = input_size - lc.offset;
    if (lc.world_rank == lc.world_size -1) chunk_len+= lc.offset;

    while (isSep(input_buffer[*read_head]) && *read_head < chunk_len) (*read_head)++;

    *word_size = 0;
    // if we are in next chunk's region
    if(*read_head == chunk_len) return;

    *word_start_index = *read_head; //start of first word; hopefully
    do{
        (*word_size)++;
        (*read_head)++;
    }while( (*word_size < lc.max_word_size) && !isSep(input_buffer[*read_head]) && (*read_head < input_size));
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
    int input_size = strlen(lc.data[0].key);

    long int chunk_len = input_size - lc.offset; // remove offset size
    if(lc.world_rank == lc.world_size -1) chunk_len+= lc.offset; // last chuck has no offset

    // skip broken word - if any
    while(!isSep(lc.data[0].key[read_head]) && read_head < chunk_len) read_head++;

    while (read_head < chunk_len) {
        find_next_word(lc.data[0].key, input_size, &read_head, & word_start_index, & word_size);
        if (word_size == 0) break;
        words[word_counter] = (char * ) malloc((word_size + 1) * sizeof(char));
        strncpy(words[word_counter], lc.data[0].key + word_start_index, word_size);
        words[word_counter][word_size] = '\0';
        word_counter++;
        if (word_counter == buffer_size) // full
            buffer_size *= 2;
            words = (char ** ) realloc(words, buffer_size * sizeof( * words));
    }

    free(lc.data[0].key);
    free(lc.data);
    lc.data = (KeyValue * ) malloc(word_counter * sizeof(KeyValue));

    #pragma omp parallel for private(i)
    for (i = 0; i < word_counter; i++) {
	lc.data[i].key = words[i];
        lc.data[i].value = 1;
    }

    lc.local_data_len = word_counter;

    // for(i = 0; i < word_counter; i ++)
    //  printf("%d: |%s| %d\n", lc.world_rank, lc.data[i].key, strlen(lc.data[i].key));

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

    if (lc.local_data_len <= 1) return; // nothing to merge

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
                lc.data[j].value+= lc.data[i].value;
                break;
            }
        }
    }
}

void convert_buckets_into_bytes(KeyValue * bucket, int bucket_size, char ** bytes, int * bytes_size) {

    int i,j;
    *bytes_size = 0;

    for (i = 0; i < bucket_size; i++) {
        *bytes_size += strlen(bucket[i].key) + 5; // 1 for space, 4 for int
    }

    if (*bytes_size == 0){
        *bytes = NULL;
        return;
    }

    *bytes = (char * ) malloc( (*bytes_size) * sizeof(char));

    j = 0;
    for (i = 0; i < bucket_size; i++) {
        int str_size = strlen(bucket[i].key);
        strncpy(*bytes + j, bucket[i].key, str_size);
        j+= str_size;
        (*bytes)[j] = '\0';
        j++;
        (*bytes)[j] = bucket[i].value & 0x000000ff;
        j++;
        (*bytes)[j] = (bucket[i].value & 0x0000ff00) >> 8;
        j++;
        (*bytes)[j] = (bucket[i].value & 0x00ff0000) >> 16;
        j++;
        (*bytes)[j] = (bucket[i].value & 0xff000000) >> 24;
        j++;
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
        int value = (recv[i] & 0x000000ff) | 
                    ((recv[i+1] & 0x000000ff) << 8) |
                    ((recv[i+2] & 0x000000ff) << 16)|
                    ((recv[i+3] & 0x000000ff) << 24);

        //printf("entry: |%s|, %d\n", buffer, value);

        // try to merge with local bucket
        int merged = 0;
        for (k = 0; k < lc.local_data_len; k++) {
            if (strlen(lc.data[k].key) == j && memcmp(lc.data[k].key, buffer, j) == 0) {
                // merge
                lc.data[k].value += value;
                // set value to -1
                recv[i] = recv[i + 1] = recv[i + 2] = recv[i + 3] = 0;

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
        int value = (recv[i] & 0x000000ff) |
                    ((recv[i+1] & 0x000000ff) << 8) |
                    ((recv[i+2] & 0x000000ff) << 16)|
                    ((recv[i+3] & 0x000000ff) << 24);

        if (value > 0) {
            char * array = (char * ) malloc(j);
            strncpy(array, buffer, j);
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
    
    int i;

    // local reduce
    reduce_local();

    // initialize bucket sizes
    int bucket_size[lc.world_size];
    for(i = 0; i < lc.world_size; i++) {
        bucket_size[i] = 0;
    }

    // find size of each bucket
    for(i = 0; i < lc.local_data_len; i++) {
        if(lc.data[i].key != NULL) {
            unsigned long word_hash = hash(lc.data[i].key);
            int receiver = word_hash % lc.world_size;
            bucket_size[receiver]++;
        }
    }

    // allocate memory for each bucket
    KeyValue ** buckets = (KeyValue ** ) malloc(lc.world_size * sizeof(KeyValue * ));
    for (i = 0; i < lc.world_size; i++) {
        buckets[i] = (KeyValue * ) malloc(bucket_size[i] * sizeof(KeyValue));
    }
    // find buckets
    int * sizes = (int * ) calloc(lc.world_size, sizeof(int));
    for (i = 0; i < lc.local_data_len; i++) {
        if(lc.data[i].key != NULL) {
            unsigned long word_hash = hash(lc.data[i].key);
            int receiver = word_hash % lc.world_size;
            buckets[receiver][sizes[receiver]++] = lc.data[i];
        }
    }

    int * sizes_bytes = (int * ) calloc(lc.world_size, sizeof(int));
    char ** send_bytes = (char ** ) malloc(lc.world_size * sizeof(char * ));

    // convert buckets to byte arrays
    for (i = 0; i < lc.world_size; i++) {
        convert_buckets_into_bytes(buckets[i], sizes[i], & send_bytes[i], & sizes_bytes[i]);
    }

    // free all unnecessary memory
    for(i = 0; i < lc.local_data_len; i++){
        if(lc.data[i].key != NULL) {
            free(lc.data[i].key);
        }
    }
    free(lc.data);
    lc.data = NULL;
    lc.local_data_len = 0;

    for (i = 0; i < lc.world_size; i++) {
        free(buckets[i]);
    }
    free(buckets);

    int * recv_sizes_bytes = (int * ) calloc(lc.world_size, sizeof(int));

    // exchange sizes between processes
    MPI_Alltoall(sizes_bytes, 1, MPI_INT, recv_sizes_bytes, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Request send_requests[lc.world_size];

    for (i = 0; i < lc.world_size; i++) {
        MPI_Isend(send_bytes[i], sizes_bytes[i], MPI_BYTE, i, 0, MPI_COMM_WORLD, &send_requests[i]);
    }

    char ** recv_bytes = (char ** ) malloc(lc.world_size * sizeof(char * ));
    MPI_Request recv_requests[lc.world_size];

    for (i = 0; i < lc.world_size; i++) {
        recv_bytes[i] = (char * ) malloc(recv_sizes_bytes[i] * sizeof(char));
        MPI_Irecv(recv_bytes[i], recv_sizes_bytes[i], MPI_BYTE, i, 0, MPI_COMM_WORLD, &recv_requests[i]);
    }

    int index_rec = 0;
    for (i = 0; i < lc.world_size; i ++) {
        // wait for any message to arrive, then merge
        MPI_Waitany(lc.world_size, recv_requests, &index_rec, MPI_STATUS_IGNORE);
        // printf("PROCESS %d RECEIVED FROM PROCESS %d\n", lc.world_rank, index_rec);
        merge(recv_bytes[index_rec], recv_sizes_bytes[index_rec]);
    }

    MPI_Waitall(lc.world_size, send_requests, MPI_STATUS_IGNORE);

    // free remaining pointers
    for (i = 0; i < lc.world_size; i++) {
        free(recv_bytes[i]);
        free(send_bytes[i]);
    }
    free(recv_bytes);
    free(recv_sizes_bytes);
    free(send_bytes);
    free(sizes_bytes);
    free(sizes);

}

void write_file() {

    int i;
    int local_size = 0;
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
    if(local_size==0) local_size++;
    // printf("Rank %d will print %d chars\n", lc.world_rank, local_size);

    // create local result
    char * result = (char * ) malloc((local_size+1) * sizeof(char));
    int j = 0;
    for (i = 0; i < lc.local_data_len; i++) {
        // get key value pair
        KeyValue kv = lc.data[i];
        j += sprintf( & result[j], "%s %d\n", kv.key, kv.value);

    }
    result[j] = '\0';

    if(local_size==1) {
        result[0] = '\n';
        result[1] = '\0';
    }
    //printf("Rank %d count:%d result: |%s|\n", lc.world_rank, local_size ,result);

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
    MPI_File_open(MPI_COMM_WORLD, "result", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, & lc.output_file);
    MPI_File_set_view(lc.output_file, 0, MPI_CHAR, result_datatype, "native", MPI_INFO_NULL);
    MPI_File_write_all(lc.output_file, result, local_size, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close( & lc.output_file);

    // clean up
    free(result);
    free(lc.data);

}
