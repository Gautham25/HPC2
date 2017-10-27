#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "lapacke.h"
#include "cblas.h"

double randomNumber(int ubound, int lbound){
    double s;
    s = ((double)rand()/(RAND_MAX))*(ubound-lbound);
    return s;
}

void assignMatVal(double *a ,int n, int ubound, int lbound){
    int i;
    for(i=0;i<n;i++){
        a[i] = randomNumber(ubound,lbound);
    }
}



int main()
{
    srand((double)time(NULL));
    int ubound = 100, lbound = 1;
    int arrayLen[] = {1000};//,2000,3000,4000,5000};
    int size = (sizeof(arrayLen)/sizeof(arrayLen[0]));
    int n,j,i;
    for(j=0;j<size;j++){
        int n = arrayLen[i];
        struct timespec tstart={0,0},tend={0,0};
        char TRANS = 'N';
        int INFO = n;
        int LDA = n;
        int LDB = n;
        int N = n;
        int NRHS = 1;
        int *IPIV = (int *)calloc(sizeof(int),n);
        double  *arrA, *arrB;
        arrA = (double *)calloc(sizeof(double),n*n);
        arrB = (double *)calloc(sizeof(double),n);
        assignMatVal(arrA,n*n,ubound,lbound);
        assignMatVal(arrB,n,ubound,lbound);
        // use new to allocate memory if you need large space
        // Here, we want to solve AX = b
        //    x1 + 2x2 + 3x3 = 1
        //    2x1 + x2 + x3  = 1
        //    x1 + x2 + x3   = 1
        // in C, you should initialize A as:
        //  A = { 1 2 3
        //        2 1 1
        //        1 1 1 }
        // IF you use this A to call LAPACK function, it gets a wrong result

        // BUT, LAPACK need the A to store in COLUMN-order
        // SO, we initial A as (for the same system):
        //  A' = { 1 2 1
        //         2 1 1
        //         3 1 1 }
        // correct solution = {0 2 -1}'
        // LU factorization

        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   a    = 1.0;
        clock_gettime(CLOCK_MONOTONIC,&tstart);
        LAPACK_dgetrf(&N,&N,arrA,&LDA,IPIV,&INFO);
        clock_gettime(CLOCK_MONOTONIC,&tend);
        double time = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        // This function solve the Ax=B directly
        //dgetrs_(&TRANS,&N,&NRHS,A,&LDA,IPIV,B,&LDB,&INFO);

        // change the order of B according to IPIV[] from LU factorization

        for(i = 0; i < N; i++)
        {
            double tmp = arrB[IPIV[i]-1];
        	arrB[IPIV[i]-1] = arrB[i];
        	arrB[i] = tmp;
        }

        // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,arrA, &N, arrB, &N);
        UPLO = 'U';
        DIAG = 'N';
        // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,arrA, &N, arrB, &N);

        printf("print the result : {\n");
        for (i=0;i<N;i++)
        {
    	       printf("%f ",arrB[i]);
        }
        printf("}\n");
        printf("Size N = %d\n",arrayLen[i]);
        printf("Time Taken = %.5f seconds\n",time);
        double gflops = (2*pow(n,3))/(3*time*pow(10,9));
        printf("\nPerformance in GFLOPS = %f\n",gflops);

    }
    return 0;
}
