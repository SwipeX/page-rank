#ifndef PAGERANK_UTIL_H__
#define PAGERANK_UTIL_H__
#include <stdlib.h>

//// general sparse matrix representation and operations////
struct SparseMatrixHandle {
	int rowdim;
	int coldim;
	int nnz;
	int * row;
	int * col;
	double * val;
};

typedef struct SparseMatrixHandle *SparseMatrix;

int makeShapeSM(SparseMatrix S, int m, int n, int nnz);
// Arrays are allocated but not filled.
// Returns 1 if success, 0 if not.

int readSM(SparseMatrix S) ;
// read matrix from standard input.
// file format is a line with n,m, and nnz,
// then nnz lines with i, j, v, where v is the (i,j) entry.
// for example: 2 2 3 / 0 0 1 / 0 1 1 / 1 0 1 /.

int writeSM(SparseMatrix S, int flag) ;
// write matrix to standard output.

//// vector operations ////

double vecdistsq(double* x, double *y, int n) ;
// For length n vectors x, y, return |x-y|^2.

double vecsum(double *x, int n) ;
// sum of entries of x.

void printvec(double *x, int n) ;

//// link matrix test case constructor 

int randomLM(SparseMatrix L, int n) ;
// Random sparse adjacency matrix. 
// Nonzeroes per row meets a power law.
// Return 1 if success, 0 if not.

/////// implementations ///////////////

int makeShapeSM(SparseMatrix S, int m, int n, int nnz) {
// Arrays are allocated but not filled.
// Returns 1 if success, 0 if not.
	S->rowdim = m;
	S->coldim = n;
	S->nnz = nnz;
	S->row = (int *) malloc(nnz*sizeof(int));
	S->col = (int *) malloc(nnz*sizeof(int));
	S->val = (double *) malloc(nnz*sizeof(double));
	if ( S->row == NULL || S->col == NULL || S->val == NULL)
		return 0;
	return 1;
}

int readSM(SparseMatrix S) {
	int m, n, nnz;
	scanf("%d %d %d\n", &m, &n, &nnz);
	makeShapeSM(S, m, n, nnz);
	int *rp = S->row;
	int *cp = S->col;
	double *vp = S->val;
	for ( ; rp < S->row + S->nnz; ++rp, ++cp, ++vp)
		scanf("%d %d %g\n", rp, cp, vp);
	return 1;
}

int writeSM(SparseMatrix S, int flag) {
// write matrix to standard output.  Whole matrix if flag = 0.  Flag at each end if flag>0.
	printf("%d %d %d\n", S->rowdim, S->coldim, S->nnz);
	int k, nnz = S->nnz;
	if (flag == 0) {
		for (k = 0; k < nnz; ++k) printf("%d %d %g\n", S->row[k], S->col[k], S->val[k]);
		return 0;
	}
	// print first flag entries
	int m = (nnz <= flag ? nnz : flag);
	for (k = 0; k < m; ++k) printf("%d %d %g\n", S->row[k], S->col[k], S->val[k]);
	// print last flag entries
	if (nnz > flag) {
		m = (nnz <= 2*flag ? flag : nnz - flag);
		for (k = m; k < nnz; ++k) printf("%d %d %g\n", S->row[k], S->col[k], S->val[k]);
	}
}

double vecdistsq(double* x, double *y, int n) {
// For length n vectors x, y, return |x-y|^2.
	int i;
	double diff, ans = 0;
	for (i = 0; i < n; ++i) { diff = x[i]-y[i]; ans += diff*diff; }
	return ans;
}

double vecsum(double *x, int n) {
// sum of entries of x.
	int i;
	double ans = 0;
	for (i = 0; i < n; ++i) { ans += x[i]; }
	return ans;
}

void printvec(double *x, int n) {
	printf("(");
	// print first 6 entries
	int i, m = (n <= 6 ? n : 6);
	for (i = 0; i < m; ++i) printf("%g, ", x[i]);
	// print last 6 entries
	if (n > 6) {
		printf("\n");
		m = (n <= 12 ? 6 : n - 6);
		for (i = m; i < n; ++i) printf("%g, ", x[i]);
	}
	printf("), sum is %g\n", vecsum(x, n));
}

void printindexvec(int *x, int n) {
	printf("(");
	int i, m = (n < 20 ? n : 20);
	for (i = 0; i < m; ++i) printf("%d, ", x[i]);
	printf(")\n");
}

int intcomp(const void *a, const void *b) 
{ int aa = *(int*)a, bb = *(int*)b; return aa < bb ? -1 : (aa == bb ? 0 : 1);}

int randomLM(SparseMatrix L, int n) {
// Random sparse adjacency matrix. 
// Nonzeroes per row meets a power law.
// Return 1 if success, 0 if not.
	int i,j,k,m;
/*
	int i,j,k, m = n, lg = 0;
	while(1 < m) { m >>= 1; lg += 1; } 
	// now 2^lg >= n > 2^(lg-1) (for n > 0).  Thus lg = ceil(log(n)).
	int nnz = n*lg*a;
	if (nnz > n*n) nnz = n*n;
	int q = 4, step = 4, m; // step is a power of 2.
	nnz = 0; 
	while (q + step < n) { nnz += step*(n/step); q += step; step *= 2; }
	// example: nnz, q,step,nnz for n = 18
	0 4 4
	16 8 8
	32 12 16
	nnz += n
--
*/
	int p2 = 1; while (p2 <= n + 1) p2 <<= 1; p2 >> 1; // 1 <= p2-1 <= n, if n > 0
	int np = n - (p2 - 1);  // p2-1 = 2^(l+1)-1
	int nnz = 0; while (p2 > 1) { p2 /= 2; nnz += p2*(n/p2); } nnz += np;

	k = nnz; // index that will step down to 0.
	makeShapeSM(L, n, n, nnz);

	// a permutation of the col indices
	int *cols = (int *) malloc(n*sizeof(int));
	for (i = 0; i < n; ++i) cols[i] = i;
		for (j = 0; j < n; ++j) { // randomize 
			int jj = random()%n; 
			int tmp = cols[j]; cols[j] = cols[jj]; cols[jj] = tmp;
		}

	//int l = n/2, 
	int q = 1; p2 = 1;
	for (i = 0; i < n; ++i) {
		if ( i >= q ) { p2 <<= 1; q += p2; if (p2 > n) p2 = n; }
		m = n/p2;
/*
		//int l = 1, key = 1 + random()%n;
		//while(key % 2 == 0) { key >>= 1; l *= 2; }
		//row i gets l*a entries;
		m = (m < 1 ? 1 : (m > n ? n : m));

		// I have k positions to fill  
		// I want to do at least 1 and at most n per row
		// I have n-i rows remaining.
		// I require n-i <= k <= n*(n-i).
		// Initially I have n-0 <= nnz-0) <= n*(n-0).
		// When n-i <= k <= n*(n-i), at the next step
		// I require n-(i+1) <= k-m <= n*(n-(i+1)).
		// This translates to k - (n-(i+1)) >= m >= k - n*(n-(i+1)) 

		int mm = k - (n-(i+1));
		if (mm < m) {// m is too big
			m = (mm > n ? n : mm); // this case should not happen: (mm < 1 ? 1 : mm));
		}
		mm = k - n*(n-(i+1));
		if (m < mm) {// m is too small
			m = (mm < 1 ? 1 : mm); // this case should not happen: (mm > n ? n : mm));
		}
*/

		int start = random()%n;
		if (start + m > n) start = n - m;
		for (j = start; j < start+m; ++j) { // randomize 
			int jj = random()%n; 
			int tmp = cols[j]; cols[j] = cols[jj]; cols[jj] = tmp;
		}
		qsort(cols+start, m, sizeof(int), intcomp);
		
		// fill the row
		for (j = start; j < start+m && k > 0; ++j, --k) {
			L->row[nnz-k] = i;
			L->col[nnz-k] = cols[j];
			L->val[nnz-k] = 1; 
		}
		//printf("row %d, m was %d, of nnz = %d, k is %d\n", i, m, nnz, k);
	}
	free(cols);
}

#endif // PAGERANK_UTIL_H__