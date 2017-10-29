#include <stdio.h>
#include <stdlib.h>
#include<time.h>
#include<math.h>
#include "lapacke.h"
#include "blas.h"

double Random_gen( int upper, int lower)
{
    double s;
    s = ((double)rand()/(RAND_MAX))*(upper-lower);
    return s;
}


//to implement transpose as fortran has column wise implementation and the program is implemented in row wise implementation
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
void mydgetrf(double *a,int *pvt,int n,int block,double *tempv){
    int maxind,temps,ib,end,i,t,k,j,p,q,l,m;
    double max,sum;
    double *ll;
    for(ib=0;ib<n;ib+=block)
    {
        end=(block-1)+ib;
        for(i=ib;i<=end;i++)
        {
           maxind=i;
           max=abs(a[i*n+i]);
           for(t=i+1;t<n;t++)
           {

                if(abs(a[t*n+i])>max)
                 {
                     maxind=t;
                     max=abs(a[t*n+i]);

                 }
           }
           if(max==0.0)
           {

            printf("LU factorization failed:coefficient matrix is singular");
            return;
           }
           else
           {

                if(maxind != i)
                {   //save pivoting information
                    temps = pvt[i];
                    pvt[i] = pvt[maxind];
                    pvt[maxind] = temps;

                    for(k=0;k<n;i++)
                    {
                        tempv[k]=a[i*n+k];
                        a[i*n+k]=a[maxind*n+k];
                        a[maxind*n+k]=tempv[k];
                    }
                }
           }

            // for(j=i+1;j<n;j++)
            // {
            //     a[j*n+i]=a[j*n+i]/a[i*n+i];
            //     for(k=i+1;k<=end;k++){
            //         a[j*n+k]=a[j*n+k]-(a[j*n+i]*a[i*n+k]);
            //     }
            // }
        }
        //factorizing

        //ll inverse
        // ll = (double*)calloc(sizeof(double), block*block);
        //     p=0;q=0;
        //     for(l=ib;l<=end;l++){
        //         for(m=ib;m<=end;m++){
        //             if(l>m){
        //                 ll[p*block+q] = a[l*n+m] * (-1);
        //             }
        //             else if(l==m){
        //                 ll[p*block+q] = 1;
        //             }
        //             else{
        //                 ll[p*block+q] = 0;
        //             }
        //             q++;
        //         }
        //         p++;
        //         q=0;
        //     }
        //     p=0;q=0;
        //     for(j=ib;j<=end;j++){
        //         for(k=end+1;k<n;k++){
        //                 sum=0.0;
        //             for(m=ib;m<=end;m++){
        //                 sum+= ll[p*block+q] * a[m*n+k];
        //                 q++;
        //             }
        //             a[j*n+k]=sum;
        //             q=0;
        //         }
        //         p++;
        //         q=0;
        //     }
        //     for(j=end+1;j<n;j++){
        //         for(k=end+1;k<n;k++){
        //             double store=0.0;
        //             for(l=ib;l<=end;l++){
        //                 store+=a[j*n+l]*a[l*n+k];
        //             }
        //             a[j*n+k]-=store;
        //         }
        //     }
        //     free(ll);
    }
}
//forward substitution
void mydtrsm_f(int n, double *a, double *B, int *pvt, double *x, double *y){
    double sum;
    int i,k;
    y[0] = B[pvt[0]];
    for(i=1;i<n;i++){
        sum = 0.0;
        for(k=0;k<i;k++){
            sum += y[k]*a[i*n+k];
        }
        y[i] = B[pvt[i]]-sum;
    }
}
//backward substitution
void mydtrsm_b(int n, double *a, double *B, int *pvt, double *x, double *y){
    double sum;
    int i,k;
    x[n-1] = y[n-1]/a[(n-1)*n+(n-1)];
    for(i=n-2;i>=0;i--){
        sum=0.0;
        for(k=i+1;k<n;k++){
            sum+= x[k]*a[i*n+k];
        }
        x[i] = (y[i]-sum)/a[i*n+i];
    }
}


int main()
{
    srand((double)time(NULL));
    int *pvt,*IPIV;
    int n,k,i,j;
    double ran=Random_gen(10,1);
    double *a,*B,*a1,*B1,*x,*y,*tempv;
    double difference,error=0.0;
    double gflops,cpu_time;
    struct timespec cstart = {0,0}, cend ={0,0};
    int block[1]={2};//50,100,200,500};
    for(n=100;n<200;n=n+1000)//loop for n values
    {
        for(k=0;k<1;k++)// loop for block_size
        {
            a=(double *) calloc(sizeof(double), n*n);
            B=(double *) calloc(sizeof(double), n);
            a1=(double *) calloc(sizeof(double), n*n);
            B1=(double *) calloc(sizeof(double), n);
            pvt=(int *) calloc(sizeof(int), n);
            y=(double *) calloc(sizeof(double), n);
            x=(double *) calloc(sizeof(double), n);
            tempv=(double *) calloc(sizeof(double), n);

            for(i=0;i<n;i++)
                for(j=0;j<n;j++)
                {
                      a[i*n+j]=Random_gen(10,1);
                      a1[i*n+j]=a[i*n+j];
                }
            for(i=0;i<n;i++){
                B[i]=Random_gen(10,1);
                B1[i]=B[i];
                pvt[i]=i;
            }

            clock_gettime(CLOCK_MONOTONIC, &cstart);
            mydgetrf(a,pvt,n,block[k],tempv);
            clock_gettime(CLOCK_MONOTONIC, &cend);
            cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
            printf("\nCPU time for LU factorization n=%d is %f",n,cpu_time);
            gflops=(2*pow(n,3))/(3*cpu_time*pow(10,9));
            printf("\nthe gflops used are=%f",gflops);
            //mydtrsm_f(n,a,B,pvt,x,y);
            //mydtrsm_b(n,a,B,pvt,x,y);
        }
        IPIV = (int *)calloc(sizeof(int),n);
        char    TRANS = 'N';
        int     NRHS = 1;
        int     INFO = n;
        int     LDA = n;
        int     LDB = n;
        //LU factorization
        transpose(a1,n);
        clock_gettime(CLOCK_MONOTONIC, &cstart);
        LAPACK_dgetrf(&n,&n,a1,&LDA,IPIV,&INFO);
        clock_gettime(CLOCK_MONOTONIC, &cend);
        cpu_time=((double)cend.tv_sec + 1.0e-9*cend.tv_nsec) - ((double)cstart.tv_sec + 1.0e-9*cstart.tv_nsec);
        printf("\nCPU time for LU factorization n=%d is %f",n,cpu_time);
        gflops=(2*pow(n,3))/(3*cpu_time*pow(10,9));
        printf("\nthe gflops used are=%f",gflops);
        char     SIDE = 'L';
        char     UPLO = 'L';
        char     DIAG = 'U';
        int      M    = 1;
        double   b    = 1.0;

        for(i = 0; i < n; i++)
        {
            double tmp = B1[IPIV[i]-1];
            B1[IPIV[i]-1] = B1[i];
            B1[i] = tmp;
        }
        // forward  L(Ux) = B => y = Ux
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&b,a1, &n, B1, &n);
        UPLO = 'U';
        DIAG = 'N';
        // backward Ux = y
        dtrsm_(&SIDE,&UPLO,&TRANS,&DIAG,&n,&M,&b,a1, &n, B1, &n);
        // for(i=0;i<n;i++){
        //    difference=(abs)(B1[i]-x[i]);
        //    if(difference>error)
        //         error=difference;
        // }
        printf("\n the error value for n=%d is %f ",n,error);

    }
    return 0;
}
