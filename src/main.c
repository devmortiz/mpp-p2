#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <sys/time.h>
#include <mpi.h>

#include "../include/io.h"

#define MODO 2 // 0 = síncrona, 1 = asíncrona, 2 = colectiva

extern double aplicar_mh(const double *, int, int, int, int, int *, int, char **, int, int, int);

static double mseconds() {
	struct timeval t;
	gettimeofday(&t, NULL);
	return t.tv_sec*1000 + t.tv_usec/1000;
}

int main(int argc, char **argv)
{
	int np, rank;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

//	Check Number of Input Args
	if(argc < 4) {
		fprintf(stderr,"Ayuda:\n"); 
		fprintf(stderr,"  ./programa n m nGen tamPob\n");
		return(EXIT_FAILURE);
	}
	
	int n = atoi(argv[1]);
	int m = atoi(argv[2]);
	int n_gen = atoi(argv[3]);
	int tam_pob = atoi(argv[4]);

	
//	Check that 'm' is less than 'n'
	assert(m < n);
	
//	Generate matrix D with distance values among elements
	double *d = read_distances(n);
	MPI_Status status;
	MPI_Request p_requests[np];

	if (MODO == 2) {
		MPI_Bcast(d, ((n*n-n)/2), MPI_DOUBLE, 0, MPI_COMM_WORLD);
		//printf("Matriz: %f, %f, %f\n", d[0], d[1], d[2]);
	}
	else {
		if(rank == 0) {
			for (int p = 1; p < np; p++) {
				if (MODO == 0)
					MPI_Send(d, ((n*n-n)/2), MPI_DOUBLE, p, 0, MPI_COMM_WORLD); // repartir el resto
				else {
					MPI_Isend(d, ((n*n-n)/2), MPI_DOUBLE, p, 0, MPI_COMM_WORLD, &p_requests[rank]); // repartir el resto
					MPI_Wait(&p_requests[rank], &status);
				}
			}

		} else {
			if (MODO == 0)
				MPI_Recv(d, ((n*n-n)/2), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &status);
			else {
				MPI_Irecv(d, ((n*n-n)/2), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, &p_requests[rank]);
				MPI_Wait(&p_requests[rank], &status);
			}
		}
	}
	
	//printf("Proceso %d, matriz: %f, %f, %f", rank, d[0], d[1], d[2]);
	
	#ifdef DEBUG
//		print_distances(d, n);
	#endif
	
//	Allocate memory for output data
	int *sol = (int *) malloc(m * sizeof(int));
	
	#ifdef TIME
		double ti = mseconds();
	#endif
	
//	Call Metaheuristic
	double value = aplicar_mh(d, n, m, n_gen, tam_pob, sol, argc, argv, np, rank, MODO);
	
	#ifdef TIME
		double tf = mseconds();
		printf("Execution Time: %.2lf sec\n", (tf - ti)/1000);
	#endif
		
	#ifdef DEBUG
		print_solution(n, m, sol, value);
	#endif

// Free Allocated Memory
  free(sol);
  free(d);

	return(EXIT_SUCCESS);
}
