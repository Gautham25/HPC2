#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include "lapacke.h"
#include "blas.h"

double randomNumber(int ubound, int lbound){
    double s;
    s = ((double)rand()/(RAND_MAX))*(ubound-lbound);
    return s;
}

void copyMatrix(double *a, double *b, int n){
    int i,j;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            b[i*n+j] = a[i*n+j];
        }
    }
}


void mydgetrf(double *arrA,int *pvt, double *tempv, int n){
    int i,t,j,k,maxind,temps;
    double max;
    for(i=0;i<n-1;i++){
        maxind = i;
        max=abs(arrA[i*n+i]);
        for(t=i+1;t<n;t++){
            if(abs(arrA[t*n+i])>max){
                maxind = t;
                max = abs(arrA[t*n+i]);
            }
        }
        if(max==0){
            printf("LU factorization failed: coefficient matrix is singular\n");
            return;
        }
        else{
            if(maxind != i){
                //Save pivoting information
                temps = pvt[i];
                pvt[i] = pvt[maxind];
                pvt[maxind] = temps;
                //Swap rows
                for(k=0;k<n;k++){
                    tempv[k] = arrA[i*n+k];
                    arrA[i*n+k] = arrA[maxind*n+k];
                    arrA[maxind*n+k] = tempv[k];
                }
            }
        }
        for(j=i+1;j<n;j++){
            arrA[j*n+i] = arrA[j*n+i]/arrA[i*n+i];
            for(k=i+1;k<n;k++){
                arrA[j*n+k] = arrA[j*n+k] - arrA[j*n+i] * arrA[i*n+k];
            }
        }
    }
}

void mydtrsm_f(int n, double *arrA, double *arrB, int *pvt, double *x, double *y){
    double sum = 0.0, temp;
    int i,k;
    printf("HELLO Begin\n");
        y[0] = arrB[pvt[0]];
        for(i=1;i<n;i++){
            sum = 0.0;
            for(k=0;k<i;k++){
                sum += y[k]*arrA[i*n+k];
            }
            y[i] = arrB[pvt[i]]-sum;
        }
        printf("HELLO\n");
}

void mydtrsm_b(int n, double *arrA, double *arrB, int *pvt, double *x, double *y){
    double sum = 0.0, temp;
    int i,k;
    x[n-1] = y[n-1]/arrA[(n-1)*n+(n-1)];
    for(i=n-2;i>=0;i--){
        sum=0.0;
        for(k=i+1;k<n;k++){
            sum+= x[k]*arrA[i*n+k];
        }
        x[i] = (y[i]-sum)/arrA[i*n+i];
    }
}

void transpose(double *a, int n){
    int i,j;
    double temp;
    for(i=0;i<n;i++){
        for(j=i;j<n;j++){
            temp = a[i*n+j];
            a[i*n+j] = a[j*n+i];
            a[j*n+i] = temp;
        }
    }
}

void assignMatVal(double *a ,int n, int ubound, int lbound){
    int i;
    for(i=0;i<n;i++){
        a[i] = randomNumber(ubound,lbound);
    }
}

double checkCorrectness(double *a, double *b, int n){
    int i,j;
    double error = 0.0;
    for(i=0;i<n;i++){
        if(error < abs(a[i]-b[i]))
            error = abs(a[i]-b[i]);
    }
    printf("Error = %f\n",error);
    printf("\n");
}

void printArray(double *a, int n, int d){
    int i,j;
    if(d==2){
        for(i=0;i<n;i++){
            for(j=0;j<n;j++){
                printf("%f ",a[i*n+j]);
            }
            printf("\n");
        }
    }
    else{
        for(i=0;i<n;i++){
            printf("%f ",a[i]);
        }
        printf("\n");
    }
}

int main()
{
    srand((double)time(NULL));
    int ubound = 100, lbound = 1;
    double random = randomNumber(ubound,lbound);
    double time,gflops;
    int arrayLen[] = {3};//,2000,3000,4000,5000};
    int size = (sizeof(arrayLen)/sizeof(arrayLen[0]));
    int n,j,i,k;
    printf("Using LAPACK Library\n");
    for(j=0;j<size;j++){
        int n = arrayLen[j];
        struct timespec tstart={0,0},tend={0,0};
        char TRANS = 'N';
        int INFO = n;
        int LDA = n;
        int LDB = n;
        int N = n;
        int NRHS = 1;
        int *IPIV = (int *)calloc(sizeof(int),n);
        double  *arrA, *arrA1, *arrB, *arrB1, *x, *y, *abk, *tempv;
        int *pvt;
        arrA = (double *)calloc(sizeof(double),n*n);
        arrA1 = (double *)calloc(sizeof(double),n*n);
        arrB = (double *)calloc(sizeof(double),n);
        arrB1 = (double *)calloc(sizeof(double), n);
        tempv = (double *)calloc(sizeof(double),n);
        assignMatVal(arrA,n*n,ubound,lbound);
        copyMatrix(arrA,arrA1,n);
        assignMatVal(arrB,n,ubound,lbound);
        for(k=0;k<n;k++){
            arrB1[k] = arrB[k];
        }
        abk = (double *)calloc(sizeof(double), n*n);
        x = (double *)calloc(sizeof(double), n);
        y = (double *)calloc(sizeof(double), n);
        pvt = (int *)calloc(sizeof(int), n);
        for(k=0;k<n;k++){
            pvt[k]=k;
        }
        transpose(arrA,n);
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
        printf("\nLAPACK LIBRARY\n");
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
        printf("\n");
        printArray(arrB,n,1);
        printf("\n");
        // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,arrA, &N, arrB, &N);

        // printf("print the result : {\n");
        // for (i=0;i<N;i++)
        // {
    	//        printf("%f ",arrB[i]);
        // }
        printf("Size N = %d\n",n);
        // printf("Time Taken = %.5f seconds\n",time);
        // double gflops = (2*pow(n,3))/(3*time*pow(10,9));
        // printf("\nPerformance in GFLOPS = %f\n",gflops);
        printf("\n");
        clock_gettime(CLOCK_MONOTONIC,&tstart);
        mydgetrf(arrA1,pvt, tempv,n);
        clock_gettime(CLOCK_MONOTONIC,&tend);
        mydtrsm_f(n,arrA1,arrB1,pvt,x,y);
        printf("\n");
        printArray(y,n,1);
        mydtrsm_b(n,arrA1,arrB1,pvt,x,y);
        // printf("\n");
        // printf("\n");
        // printf("MYDGETRF VERSION\n");
        // time = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        // gflops = (2*pow(n,3))/(3*time*pow(10,9));
        // printf("Time Taken = %.5f seconds\n",time);
        // printf("\nPerformance in GFLOPS = %f\n",gflops);
        printf("\n");
        //printArray(IPIV,n,1);
        printf("\n");
        //printArray(x,n,1);
        free(arrA);
        free(arrB);
        free(arrA1);
        free(arrB1);
        free(pvt);
        free(x);
        free(y);
        free(abk);
        free(IPIV);
        free(tempv);
    }
    return 0;
}
