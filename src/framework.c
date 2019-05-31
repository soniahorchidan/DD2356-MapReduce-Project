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
    char* filename;
    MPI_File input_file, output_file;

    // Datatype for reading/writing
    MPI_Datatype read_type;
    MPI_Datatype write_type;

    // local data
    KeyValue * big_bucket;
    KeyValue * small_bucket;

    // ranks and sizes
    int world_rank, world_size;

    // input size
    long input_file_len;
    //local sizes
    int big_bucket_size;
    int small_bucket_size;
    // max word size
    int max_word_size;
    // chunk size per iteration
    int chunk_size;
    // offset in the end of chunk
    int chunk_offset;
    // iteration counter
    int iteration_counter;
    // number of iterations until finish
    int total_iterations;
    // number of rows in the input matrix
    int rows;
    // whethe the process has read data
    int offset; // to be removed, kept for compiling
}
LocalConfig;

LocalConfig lc = {.chunk_size = 67708864, .max_word_size = 16}; //67108864

int  read_file(char * path_to_file) {

    lc.filename = path_to_file;

    // get world details
    MPI_Comm_rank(MPI_COMM_WORLD, & lc.world_rank);
    MPI_Comm_size(MPI_COMM_WORLD, & lc.world_size);

    // read input file size
    if (lc.world_rank == 0) {
        MPI_File_open(MPI_COMM_SELF, lc.filename, MPI_MODE_RDONLY, MPI_INFO_NULL, & lc.input_file);
        MPI_Offset size;
        MPI_File_get_size(lc.input_file, & size);
        MPI_File_close( & lc.input_file);
        lc.input_file_len = size - 1; // exclude EOF
        printf("Input file contains %ld chars.\n", lc.input_file_len);
    }

    // broadcast input size
    MPI_Bcast( & lc.input_file_len, 1, MPI_LONG, 0, MPI_COMM_WORLD);

    lc.big_bucket_size = 0;
    lc.iteration_counter = 0;
    lc.rows = lc.input_file_len / lc.chunk_size;
    lc.rows += (lc.input_file_len % lc.chunk_size== 0)? 0 : 1;
    lc.total_iterations = lc.rows / lc.world_size;
    lc.total_iterations += (lc.rows % lc.world_size == 0)? 0 : 1;

    MPI_File_open(MPI_COMM_WORLD, lc.filename, MPI_MODE_RDONLY, MPI_INFO_NULL, & lc.input_file);

    if(lc.world_rank == 0) {
        printf("Number of iterations: %d.\n", lc.total_iterations);
        printf("Number of rows: %d.\n", lc.rows);
    }

    return lc.total_iterations;
}

void read_chunk() {

    // calculate local data size
    int read_len = lc.chunk_size + lc.max_word_size;
    // first round, first process
    int frfp = (lc.iteration_counter == 0) && (lc.world_rank == 0);
    int last_chunk_read_by = (lc.rows % lc.world_size) -1;
    if(last_chunk_read_by == -1) last_chunk_read_by = lc.world_size -1;
    // last round, last process (this is not world_size -1 )
    int lrlp = (lc.iteration_counter== lc.total_iterations-1) && (lc.world_rank == last_chunk_read_by);
    // last round, second last process
    int lrslp = (lc.iteration_counter == lc.total_iterations-1) && (lc.world_rank == last_chunk_read_by -1);
    // if proc should read at last round
    int should_not_read = (lc.iteration_counter== lc.total_iterations-1) && (last_chunk_read_by < lc.world_rank);

    // first process in first round cannot read one byte from prev chunk
    if(frfp) {
        read_len -=1;
    }

    // last process in the last round
    if(lrlp){
        read_len = 1 + (lc.input_file_len % lc.chunk_size);
        if(read_len == 1) read_len+= lc.chunk_size;
    }

    // the second last process in the last round should read the remaining bytes
    if(lrslp){
        int front = lc.max_word_size -1;
        int rem = lc.input_file_len % lc.chunk_size;
        read_len -= front;
        read_len += (((front)<(rem))?(front):(rem));
    }

    // save chunk front offset
    lc.chunk_offset = read_len - lc.chunk_size -1;
    if(frfp) lc.chunk_offset++;
    if(lc.chunk_offset < 0) lc.chunk_offset = 0;

    MPI_Offset file_offset = (lc.world_rank-1)+(lc.iteration_counter*lc.world_size);
    file_offset *= (long)lc.chunk_size;
    file_offset += lc.chunk_size -1;

    if(frfp){
        file_offset = 0;
    }

    if(should_not_read) { //there is nothing to read
        // dummy values
        file_offset = 0;
        read_len = 1;
    }

    // alloc space for reading
    char *p, *p_dummy;

    if(frfp){
        p_dummy = (char * ) malloc((read_len+2) * sizeof(char));
        p = &p_dummy[1]; // we save the first slot for a separator
        p_dummy[0] = ' ';
        p_dummy[read_len+1] = '\0';
    }else{
        p = (char * ) malloc((read_len+1) * sizeof(char));
        p[read_len] = '\0';
    }

    // read file
    MPI_File_read_at_all(lc.input_file, file_offset, p, read_len, MPI_CHAR, MPI_STATUS_IGNORE);

    if(frfp) {
        p = p_dummy;
    }

    lc.iteration_counter++;
    if(lc.iteration_counter == lc.total_iterations) {
        MPI_File_close(&lc.input_file);
    }

    //printf("Iter:%d Rank:%d read: |%s|\n",lc.iteration_counter,lc.world_rank ,p);

    if(should_not_read){
        free(p);
        lc.small_bucket = NULL;
        lc.small_bucket_size = 0;
        return;
    }
    lc.small_bucket = (KeyValue * ) malloc(sizeof(KeyValue));
    lc.small_bucket[0].key = p;
    lc.small_bucket[0].value = 1;
    lc.small_bucket_size = 1;
}


int isSep(char c) {
    return !isdigit(c) && !isalpha(c);
}

void find_next_word(char * input_buffer, int input_size, int * read_head, int * word_start_index, int * word_size) {
    int chunk_len = input_size - lc.chunk_offset;

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

    if (lc.small_bucket_size == 0) return;
    int buffer_size = 100;
    int word_start_index;
    int word_size;
    char ** words = (char ** ) malloc(buffer_size * sizeof(char * ));
    int i;
    int word_counter = 0;
    int read_head = 0;
    int input_size = strlen(lc.small_bucket[0].key);
    int chunk_len = input_size - lc.chunk_offset; // remove offset size

    // skip broken word - if any
    while(!isSep(lc.small_bucket[0].key[read_head]) && read_head < chunk_len) read_head++;

    while (read_head < chunk_len) {
        find_next_word(lc.small_bucket[0].key, input_size, &read_head, & word_start_index, & word_size);
        if (word_size == 0) break;
        words[word_counter] = (char * ) malloc((word_size + 1) * sizeof(char));
        strncpy(words[word_counter], lc.small_bucket[0].key + word_start_index, word_size);
        words[word_counter][word_size] = '\0';
        word_counter++;
        if (word_counter == buffer_size) // full
            buffer_size *= 2;
        words = (char ** ) realloc(words, buffer_size * sizeof( * words));
    }

    free(lc.small_bucket[0].key);
    free(lc.small_bucket);
    lc.small_bucket = (KeyValue * ) malloc(word_counter * sizeof(KeyValue));

    #pragma omp parallel for private(i)
    for (i = 0; i < word_counter; i++) {
	lc.small_bucket[i].key = words[i];
        lc.small_bucket[i].value = 1;
    }

    lc.small_bucket_size = word_counter;

    //for(i = 0; i < word_counter; i ++)
    //    printf("%d: |%s| %d\n", lc.world_rank, lc.small_bucket[i].key, strlen(lc.small_bucket[i].key));

    free(words);
}

unsigned long hash(char * str) {
    unsigned long hash = 5381;
    int c;

    while ((c = * str++))
        hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

    return hash;
}


void local_reduce(KeyValue *bucket, int size) {
    // reduce sorted bucket
    
    int i = 0;          // iterate the bucket

    while(i < size) {
        int j = i + 1;
        int current_word_size = strlen(bucket[i].key);
        while(j < size && current_word_size == strlen(bucket[j].key)   
            && memcmp(bucket[i].key, bucket[j].key, current_word_size) == 0) {
            free(bucket[j].key);
            bucket[j].key = NULL;
            bucket[i].value += bucket[j].value;
            j ++;
        }
        i = j;
    }
}


void convert_buckets_into_bytes(KeyValue * bucket, int bucket_size, char ** bytes, int * bytes_size) {

    int i,j;
    *bytes_size = 4; // array starts with an int showing the number or pairs

    for (i = 0; i < bucket_size; i++) {
        *bytes_size += strlen(bucket[i].key) + 5; // 1 for space, 4 for int
    }

    *bytes = (char * ) malloc( (*bytes_size) * sizeof(char));
    j = 0;

    // append number or pairs
    (*bytes)[j]   =  bucket_size & 0x000000ff;
    (*bytes)[j+1] = (bucket_size & 0x0000ff00) >> 8;
    (*bytes)[j+2] = (bucket_size & 0x00ff0000) >> 16;
    (*bytes)[j+3] = (bucket_size & 0xff000000) >> 24;
    j+=4;


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


int compare (const void * a, const void * b) {
  KeyValue *kv1 = (KeyValue *)a;
  KeyValue *kv2 = (KeyValue *)b;
  return strcmp(kv1->key, kv2->key);
}

void sort_bucket(KeyValue *bucket, long len) {
    qsort(bucket, len, sizeof(KeyValue), compare);
}


int find_next_kv_pair(char *recv, int *current, char buffer[], int *value) {
    // takes the byte array *recv and a current index and returns the next position
    // of the index, the key of the new kv pair (in buffer variable) and the 
    // value of the new kv pair (in the value variable)
    
    int k = 0;
    while (recv[*current] != '\0') {
        buffer[k] = recv[*current];
        (*current) ++;
        k++;
    }
    (*current)++; // for '\0'
    buffer[k] = '\0';

    // read value
    *value = (recv[*current] & 0x000000ff) | 
                ((recv[(*current)+1] & 0x000000ff) << 8) |
                ((recv[(*current)+2] & 0x000000ff) << 16)|
                 ((recv[(*current)+3] & 0x000000ff) << 24);
    (*current) += 4;
    return k;
}

void merge_bucket_with_byte_array(char *recv, int num_recv) {

    int size = 0;
    int i = 0;
    int j = 4;
    int pairs_recv = (recv[0] & 0x000000ff)        |
                     ((recv[1] & 0x000000ff) << 8 )|
                     ((recv[2] & 0x000000ff) << 16)|
                     ((recv[3] & 0x000000ff) << 24);

    KeyValue *merged = (KeyValue *)malloc((lc.big_bucket_size + pairs_recv) * sizeof(KeyValue));

    char buffer[17];

    while(i < lc.big_bucket_size && j < num_recv) {
        int current = j;
        int value;
        int key_size = find_next_kv_pair(recv, &current, buffer, &value);
        int cmp = strcmp(lc.big_bucket[i].key, buffer);

        if(cmp < 0)
            merged[size ++] = lc.big_bucket[i ++];
        else if(cmp == 0) {
            merged[size] = lc.big_bucket[i ++];
            merged[size ++].value += value;
            j = current;
        } else {
            merged[size].key = (char *) malloc(key_size + 1);
            memcpy(merged[size].key, buffer, key_size);
            merged[size].key[key_size] = '\0';
            merged[size ++].value = value;
            j = current;
        }
    }

    while(i < lc.big_bucket_size)
        merged[size ++] = lc.big_bucket[i ++];

    free(lc.big_bucket);
    while(j < num_recv){
        int value;
        int key_size = find_next_kv_pair(recv, &j, buffer, &value);
        merged[size].key = (char *) malloc(key_size + 1);
        memcpy(merged[size].key, buffer, key_size);
        merged[size].key[key_size] = '\0';
        merged[size ++].value = value;
    }

    free(recv);
    lc.big_bucket_size = size;
    lc.big_bucket = (KeyValue *)realloc(merged, size * sizeof(KeyValue));
}


void reduce() {
    int i;
    // local reduce
    sort_bucket(lc.small_bucket, lc.small_bucket_size);
    local_reduce(lc.small_bucket, lc.small_bucket_size);

    // initialize bucket sizes
    int bucket_size[lc.world_size];
    for(i = 0; i < lc.world_size; i++){
        bucket_size[i] = 0;
    }

    // find size of each bucket
    for(i = 0; i < lc.small_bucket_size; i++) {
        if(lc.small_bucket[i].key != NULL) {
            unsigned long word_hash = hash(lc.small_bucket[i].key);
            int receiver = word_hash % lc.world_size;
            bucket_size[receiver]++;
        }
    }

    // allocate memory for each bucket
    KeyValue *buckets[lc.world_size];
    for (i = 0; i < lc.world_size; i++) {
        buckets[i] = (KeyValue * ) malloc(bucket_size[i] * sizeof(KeyValue));
    }
    // find buckets
    int sizes[lc.world_size];
    for(i = 0; i < lc.world_size; i++){
        sizes[i] = 0;
    }
    for (i = 0; i < lc.small_bucket_size; i++) {
        if(lc.small_bucket[i].key != NULL) {
            unsigned long word_hash = hash(lc.small_bucket[i].key);
            int receiver = word_hash % lc.world_size;
            buckets[receiver][sizes[receiver]++] = lc.small_bucket[i];
        }
    }

    int sizes_bytes[lc.world_size];
    char *send_bytes[lc.world_size];

    // convert buckets to byte arrays
    for (i = 0; i < lc.world_size; i++) {
        convert_buckets_into_bytes(buckets[i], sizes[i], & send_bytes[i], & sizes_bytes[i]);
    }

    // free all unnecessary memory
    for(i = 0; i < lc.small_bucket_size; i++){
        if(lc.small_bucket[i].key != NULL) {
            free(lc.small_bucket[i].key);
        }
    }
    free(lc.small_bucket);
    lc.small_bucket = NULL;
    lc.small_bucket_size = 0;

    for (i = 0; i < lc.world_size; i++) {
        free(buckets[i]);
    }

    int recv_sizes_bytes[lc.world_size];
    // exchange sizes between processes
    MPI_Alltoall(sizes_bytes, 1, MPI_INT, recv_sizes_bytes, 1, MPI_INT, MPI_COMM_WORLD);

    MPI_Request send_requests[lc.world_size];

    for (i = 0; i < lc.world_size; i++) {
        MPI_Isend(send_bytes[i], sizes_bytes[i], MPI_BYTE, i, 0, MPI_COMM_WORLD, &send_requests[i]);
    }

    char *recv_bytes[lc.world_size];
    MPI_Request recv_requests[lc.world_size];

    for (i = 0; i < lc.world_size; i++) {
        recv_bytes[i] = (char * ) malloc(recv_sizes_bytes[i] * sizeof(char));
        MPI_Irecv(recv_bytes[i], recv_sizes_bytes[i], MPI_BYTE, i, 0, MPI_COMM_WORLD, &recv_requests[i]);
    }

    //wait for all the small buckets to arrive
    int index_rec;
    for (i = 0; i < lc.world_size; i ++) {
        // wait for any message to arrive, then merge
        MPI_Waitany(lc.world_size, recv_requests, &index_rec, MPI_STATUS_IGNORE);
        // printf("PROCESS %d RECEIVED FROM PROCESS %d\n", lc.world_rank, index_rec);
        merge_bucket_with_byte_array(recv_bytes[index_rec], recv_sizes_bytes[index_rec]);
    }

    MPI_Waitall(lc.world_size, send_requests, MPI_STATUS_IGNORE);

    // free remaining pointers
    for (i = 0; i < lc.world_size; i++) {
        free(send_bytes[i]);
    }
}

void write_file() {

    int i;
    int local_size = 0;
    for (i = 0; i < lc.big_bucket_size; i++) {
        KeyValue kv = lc.big_bucket[i];
        // size of string
        local_size += strlen(kv.key);
        // size of count
        local_size += floor(log10(abs(kv.value))) + 1;
        // space and newline
        local_size += 2;
    }
    if(local_size==0) local_size++;
    // printf("Rank %d will print %d chars\n", lc.world_rank, local_size);

    // create local result
    char * result = (char * ) malloc((local_size+1) * sizeof(char));
    int j = 0;
    for (i = 0; i < lc.big_bucket_size; i++) {
        // get key value pair
        KeyValue kv = lc.big_bucket[i];
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
    int start_from[1] = {
        proc_offset
    };
    int local_array[1] = {
        local_size
    };
    int result_array[1] = {
        total_size
    };
    MPI_Type_create_subarray(1, result_array, local_array, start_from, MPI_ORDER_C, MPI_CHAR, & lc.write_type);
    MPI_Type_commit( & lc.write_type);
    MPI_File_open(MPI_COMM_WORLD, "result", MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, & lc.output_file);
    MPI_File_set_view(lc.output_file, 0, MPI_CHAR, lc.write_type, "native", MPI_INFO_NULL);
    MPI_File_write_all(lc.output_file, result, local_size, MPI_CHAR, MPI_STATUS_IGNORE);
    MPI_File_close( & lc.output_file);

}


