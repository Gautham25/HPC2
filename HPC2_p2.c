#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

double randomNumber(int ubound, int lbound){
    double s;
    s = ((double)rand()/(RAND_MAX))*(ubound-lbound);
    return s;
}

double checkCorrectness(double *a, double *b, int n){
    int i,j;
    double error = 0.0;
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
            if(error < abs(a[i*n+j]-b[i*n+j]))
                error = abs(a[i*n+j]-b[i*n+j]);
        }
    }
    printf("Error = %f\n",error);
}

void assignMatVal(double *a, int n, int ubound, int lbound){
    int i;
    for(i=0;i<n;i++){
        a[i] = randomNumber(ubound, lbound);
    }
}

void matInverse(double *ll, int n){
    int i,j;
    for(i=0; i<n;i++){
        for(j=0;j<n;j++){
            if(j>i){
                ll[i*n+j] = 0;
            }
            else if(j<i){
                ll[i*n+j] *= (-1);
            }
        }
    }
}

void mydgetrf(double *arrA,int *pvt, int n, int b){
    int i,t,l,m,p,q,j,k,maxind,temps,end;
    double tempv,max;
    for(ib=0;ib<n;ib+=b){
        end = ib + b-1;
        for(i=0;i<n-1;i++){
            maxind = i;
            max=abs(arrA[i*n+i]);
            for(t=i+1;t<n;t++){
                if(abs(arrA[t*n+i])>max){
                    maxind = t;
                    max = abs(arrA(t*n+i));
                }
            }
            if(max==0){
                printf("LU factorization failed: coefficient matrix is singular\n");
                return;
            }
            else{
                if(maxind != 1){
                    //Save pivoting information
                    temps = pvt[i];
                    pvt[i] = pvt[maxind];
                    pvt[maxind] = temps;
                    //Swap rows
                    for(k=i;k<n;k++){
                        tempv = arrA[i*n+k];
                        arrA[i*n+k] = arrA[maxind*n+k];
                        arrA[maxind*n+k] = tempv;
                    }
                }
            }
            double *ll = (double*)calloc(sizeof(double), b*b);
            p=0;q=0;
            for(l=ib;l<=end;l++){
                for(m=ib;m<=end;m++){
                    ll[p*b+q] = arrA[p*b+q];
                    q++;
                }
                p++;
                q=0;
            }
            matInverse(ll,b);
            p=0;q=0;
            for(j=ib;j<=end;j++){
                //arrA[j*n+i] = arrA[j*n+i]/arrA[i*n+i];
                for(k=end+1;k<n;k++){
                    arrA[j*n+k] = ll[p*b+q] * arrA[i*n+k]
                    q++;
                }
                p++;
                q=0;
            }
            for(j=end+1;j<n;j++){
                for(k=end+1;k<n;k++){
                    for(l=ib;l<=end;l++){

                    }
                }
            }
        }
    }

}

void mydtrsm(int n, double *arrA, double *arrB, int *pvt, double *z, double *y, int label){
    double sum = 0.0, temp;
    int i,k;
    if(label == 0){
        y[0] = arrB[pvt[0]];
        for(i=1;i<n;i++){
            for(k=0;k<i-1;k++){
                sum+=y[k]*arrA[i*n+k];
                y[i] = b[pvt[i]]-sum;
            }
        }
    }
    else{
        x[n-1] = y[n-1]/a[n*n+n];
        for(i=n-1;i>=0;i--){
            for(k=i+);k<n;k++){
                sum+= x[k]*a[i*n+k];
                temp = y[i]-sum;
                x[i] = temp/a[i*n+i];
            }
        }
    }
}

int main(){
    int *pvt,n,i,k;
    int ubound = 100, lbound = 0;
    int arrN = {1000,2000,3000,4000,5000};
    double *arrA, *arrB, *abk, *x, *y;
    double gflops, time;
    struct timespec tstart={0,0}, tend={0,0};
    int len = sizeof(arrN)/sizeof(arrN[0]);
    for(i=0;i<len;i++){
        n = arrN[i];
        arrA = (double *)calloc(sizeof(double), n*n);
        arrB = (double *)calloc(sizeof(double), n);
        abk = (double *)calloc(sizeof(double), n*n);
        x = (double *)calloc(sizeof(double), n*n);
        y = (double *)calloc(sizeof(double), n*n);
        pvt = (int *)calloc(sizeof(int), n);
        for(k=0;k<n;k++){
            pvt[k]=k;
        }
        assignMatVal(arrA, n*n, ubound, lbound);
        assignMatVal(arrB, n, ubound, lbound);
        clock_gettime(CLOCK_MONOTONIC, &tstart);
        mydgetrf(arrA,pvt,n);
        clock_gettime(CLOCK_MONOTONIC, &tend);
        time = ((double)tend.tv_sec + 1.0e-9*tend.tv_nsec) - ((double)tstart.tv_sec + 1.0e-9*tstart.tv_nsec);
        printf("\nCPU time for LU factorisation is %.5f\n",time);
        gflops = (2*pow(n,3))/(3*time*pow(10,9));
        printf("\nPerformance in GFLOPS = %f\n",gflops);
        mydtrsm(n,arrA,arrB,pvt,x,y,0);
        mydtrsm(n,arrA,arrB,pvt,x,y,1);
        free(arrA);
        free(arrB);
        free(pvt);
        free(x);
        free(y);
        free(abk);
    }
    return 0;
}
