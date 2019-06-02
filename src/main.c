#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

#include "framework.h"

void print_usage();

int main(int argc, char *argv[])
{
	int opt;
	int world_rank;
	int repeat = 1;

	double avg_runtime = 0.0, prev_avg_runtime = 0.0, stddev_runtime = 0.0;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	while ((opt = getopt(argc, argv, "cfr:")) != -1) {
		switch (opt) {
			case 'r':
				repeat = atoi(optarg);
				break;
			default:
				if (world_rank == 0) print_usage(argv[0]);
				MPI_Finalize();
				exit(1);
		}
	}

	if (argv[optind] == NULL) {
		if (world_rank == 0) print_usage(argv[0]);
		MPI_Finalize();
		exit(1);
	}

	int i,j;
	double start,end;
	for (i = 0; i < repeat; i++){

		MPI_Barrier(MPI_COMM_WORLD);
		start = MPI_Wtime();

		// READ
		int rounds = read_file(argv[optind]);

		for(j = 0; j < rounds; j++){

			read_chunk();		
			// MAP
			flat_map();

			// REDUCE
			reduce();
		}
		
		// WRITE
		write_file();

		MPI_Barrier(MPI_COMM_WORLD);
		end = MPI_Wtime();
		
		if (world_rank == 0) {
			printf("run %d: %f s\n", i, end - start);
		}
		prev_avg_runtime = avg_runtime;
		avg_runtime = avg_runtime + ( (end - start) - avg_runtime ) / (i + 1);
		stddev_runtime = stddev_runtime + ( (end - start) - avg_runtime) * ( (end - start) - prev_avg_runtime);
	}

	if (world_rank == 0) {
		stddev_runtime = sqrt(stddev_runtime / (repeat - 1));
		printf("duration\t= %f±%f\n", avg_runtime, stddev_runtime);
	}

	MPI_Finalize();

	return 0;
}

void print_usage(char *program) {
	fprintf(stderr, "Usage: %s -r [Repeat run number] [Input file]\n", program);
}
