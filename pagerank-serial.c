#include <stdio.h>
#include "pagerank-util.h"
#include <mpi.h>

// pagerank-straw provides straw_mvpSM(), scale(), solve().
// some subset of these functions (definitely including straw_mvpSM) will be modified for performance.
#include "pagerank-straw.h"

/////////////////////////////////////////////////////////////
// The imagined pagerank-serial.h contains my (student) serial_mvpSM and possibly other reworked functions.
//#include "pagerank-serial.h"
// ...or I put it right here:
void serial_mvpSM(double * w, SparseMatrix S, double * z, int size, int rank) {
	// this is a stub.  ..replace this with an improved serial matrix vector product.
	strawman_mvpSM(w, S, z, size, rank);
}
/////////////////////////////////////////////////////////////

void uniform(double *x, int n) {
// Vector x has length n, is already allocated.  All entries are set to 1/n.
	int i;
	double ninv = 1.0/n;
	for (i = 0; i < n; ++i) x[i] = ninv;
}

// main tests strawman version against my new serial version.
int main(int argc, char* argv[]){
	// hello
	int n = 100;
	double d = 0.85; // damping variable
	double eps = 0.000000001;
	if (argc <= 1 || argc > 4) {
		printf("usage: %s num-pages damping-factor(0.85) epsilon(0.0001)\n", argv[0]);
		return 0;
	}
	if (argc > 1) n = atoi(argv[1]);
	if (argc > 2) d = atof(argv[2]);
	if (argc > 3) eps = atof(argv[3]);

	// build vectors
	double *y0  = (double *) malloc(n*sizeof(double)); // holds initial page probabilities
	double *y  = (double *) malloc(n*sizeof(double)); // holds intermediate and ultimate probs.
	uniform(y0, n);
	printf("initial probabilities\n");
	printvec(y0, n);

	// build matrix
	struct SparseMatrixHandle SH;
	SparseMatrix S = &SH;
	randomLM(S, n); // The link matrix
	if (n <= 100) writeSM(S, 50);
	printf("dimension is %d, nnz is %d, damper is %f, epsilon is %g\n\n", S->coldim, S->nnz, d, eps);

	double *C = (double *) malloc(S->coldim*sizeof(double));
	scale(C, S); // convert link matrix to stochastic matrix (col sums are 1).

	MPI_Init(0,0);
	double start, elapsed, elapsed2;
	int iters; 
	// compute page rank
	int size;
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	
	start = MPI_Wtime();
	iters = solve(S, C, d, y0, y, eps, strawman_mvpSM, size, rank);
	elapsed = MPI_Wtime() - start;
	printf("final (page rank) probabilities, %d iterations in time %f\n", iters, elapsed);
	printvec(y, n); printf("\n");

	// good bye
	free(C); free(y); free(y0); free(SH.row); free(SH.col); free(SH.val);
	return 0;
}
