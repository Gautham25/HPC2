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

void assignMatVal(double *a, int n, int ubound, int lbound){
    int i;
    for(i=0;i<n;i++){
        a[i] = randomNumber(ubound, lbound);
    }
}

void copyMatrix(double *a, double *b, int n){
    int i,j;
    for(i=0;i<n;i++){
        b[i] = a[i];
    }
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

void mydgetrfBlock(double *arrA,int *pvt, double *tempv, int n, int b){
    int i,t,l,m,p,q,j,k,maxind,temps,end,ib;
    double max, matsum;
    double *ll;
    for(ib=0;ib<n;ib+=b){
        end = ib + b-1;
        for(i=ib;i<=end;i++){
            maxind = i;
            max=abs(arrA[i*n+i]);
            for(t=i+1;t<n;t++){
                if(abs(arrA[t*n+i])>max){
                    maxind = t;
                    max = abs(arrA[t*n+i]);
                }
            }
            if(max==0.0){
                //printf("LU factorization failed: coefficient matrix is singular\n");
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
                for(k=i+1;k<=end;k++){
                    arrA[j*n+k] = arrA[j*n+k] - arrA[j*n+i] * arrA[i*n+k];
                }
            }
        }

        ll = (double*)calloc(sizeof(double), b*b);
        //double *lln = (double*)calloc(sizeof(double), b*b);
        p=0;q=0;
        for(l=ib;l<=end;l++){
            for(m=ib;m<=end;m++){
                if(l>m){
                    ll[p*b+q] = arrA[l*n+m] * (-1);
                }
                else if(l==m){
                    ll[p*b+q] = 1;
                }
                else{
                    ll[p*b+q] = 0;
                }
                q++;
            }
            p++;
            q=0;
        }
        p=0;q=0;
        for(j=ib;j<=end;j++){
            for(k=end+1;k<n;k++){
                matsum = 0.0;
                for(m=ib;m<=end;m++){
                    matsum += ll[p*b+q] * arrA[m*n+k];
                    q++;
                }
                arrA[j*n+k] = matsum;
                q=0;
            }
            p++;
            q=0;
        }
        for(j=end+1;j<n;j++){
            for(k=end+1;k<n;k++){
                double gmat = 0.0;
                for(l=ib;l<=end;l++){
                    gmat += arrA[j*n+l] * arrA[l*n+k];
                }
                arrA[j*n+k] -= gmat;
            }
        }
        free(ll);
    }

            //free(lln);
}

void mydtrsm_f(int n, double *arrA, double *arrB, int *pvt, double *x, double *y){
    double sum = 0.0, temp;
    int i,k;
    y[0] = arrB[pvt[0]];
    for(i=1;i<n;i++){
        sum = 0.0;
        for(k=0;k<i;k++){
            sum += y[k]*arrA[i*n+k];
        }
        y[i] = arrB[pvt[i]]-sum;
    }
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

int main(){
    srand((double)time(NULL));
    int *pvt,n,i,k,m,z;
    int ubound = 100, lbound = 0;
    int arrN[] = {1000,2000,3000,4000,5000};
    int block[] = {200};
    double random = randomNumber(ubound,lbound);
    double time,gflops;
    int len = sizeof(arrN)/sizeof(arrN[0]);
    int blockLen = sizeof(block)/sizeof(block[0]);
    for(i=0;i<5;i++){
        n = arrN[i];
        struct timespec tstart={0,0},tend={0,0};
        char TRANS = 'N';
        int INFO = n;
        int LDA = n;
        int LDB = n;
        int N = n;
        int NRHS = 1;
        int *IPIV = (int *)calloc(sizeof(int),n);
        double  *arrA, *arrA1, *arrB, *arrB1, *arrA2, *arrB2, *x, *y, *abk, *tempv;
        int *pvt;
        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   a    = 1.0;
        arrA = (double *)calloc(sizeof(double),n*n);
        arrA1 = (double *)calloc(sizeof(double),n*n);
        arrB = (double *)calloc(sizeof(double),n);
        arrB1 = (double *)calloc(sizeof(double), n);

        assignMatVal(arrA,n*n,ubound,lbound);
        copyMatrix(arrA,arrA1,n*n);
        assignMatVal(arrB,n,ubound,lbound);
        // printArray(arrB,n,1);
        for(k=0;k<n;k++){
            arrB1[k] = arrB[k];
        }
        transpose(arrA,n);
        printf("\nLAPACK LIBRARY\n");
        clock_gettime(CLOCK_MONOTONIC,&tstart);
        LAPACK_dgetrf(&N,&N,arrA,&LDA,IPIV,&INFO);
        clock_gettime(CLOCK_MONOTONIC,&tend);
        time = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        // This function solve the Ax=B directly
        //dgetrs_(&TRANS,&N,&NRHS,A,&LDA,IPIV,B,&LDB,&INFO);

        // change the order of B according to IPIV[] from LU factorization

        for(z = 0; z < N; i++)
        {
            double tmp = arrB[IPIV[z]-1];
        	arrB[IPIV[z]-1] = arrB[z];
        	arrB[z] = tmp;
        }

        // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,arrA, &N, arrB, &N);
        UPLO = 'U';
        DIAG = 'N';

        // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&N,&M,&a,arrA, &N, arrB, &N);

        // printf("print the result : {\n");
        // for (i=0;i<N;i++)
        // {
    	//        printf("%f ",arrB[i]);
        // }
        printf("Size N = %d\n",n);
        printf("Time Taken = %.5f seconds\n",time);
        gflops = (2*pow(n,3))/(3*time*pow(10,9));
        printf("\nPerformance in GFLOPS = %f\n",gflops);
        // printf("\n");
        // printArray(x,n,1);
        // printf("\n");
        printf("\nBLOCKED GEPP \n");
        printf("\nSize N = %d\n",n);
        arrA2 = (double *)calloc(sizeof(double),n*n);
        arrB2 = (double *)calloc(sizeof(double),n);
        tempv = (double *)calloc(sizeof(double),n);
        x = (double *)calloc(sizeof(double), n);
        y = (double *)calloc(sizeof(double), n);
        pvt = (int *)calloc(sizeof(int), n);
        for(k=0;k<1;k++){

            for(m=0;m<n;m++){
                pvt[m]=m;
            }
            copyMatrix(arrA1,arrA2,n*n);
            copyMatrix(arrB1,arrB2,n);
            // printArray(arrB2,n,1);
            clock_gettime(CLOCK_MONOTONIC, &tstart);
            mydgetrfBlock(arrA2,pvt,tempv,n,block[k]);
            clock_gettime(CLOCK_MONOTONIC, &tend);
            printf("Block Size = %d",block[k]);
            time = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
            printf("\nCPU time for LU factorisation is %.5f\n",time);
            gflops = (2*pow(n,3))/(3*time*pow(10,9));
            printf("\nPerformance in GFLOPS = %f\n",gflops);
            mydtrsm_f(n,arrA2,arrB2,pvt,x,y);
            mydtrsm_b(n,arrA2,arrB2,pvt,x,y);
            // printf("\n");
            // printArray(x,n,1);
            // printf("\n");
            checkCorrectness(arrB,x,n);

        }

        printf("\n");
    }
    // free(pvt);
    // free(x);
    // free(y);
    // free(arrA2);
    // free(arrB2);
    // free(tempv);
    // free(arrA);
    // free(arrB);
    // free(arrA1);
    // free(arrB1);
    // free(IPIV);
    return 0;
}
