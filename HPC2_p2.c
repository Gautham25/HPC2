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

void getCofactor(double* A, double* temp, int p, int q, int n)
{
    int i = 0, j = 0;
    int row, col;

    // Looping for each element of the matrix
    for (row = 0; row < n; row++)
    {
        for (col = 0; col < n; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != p && col != q)
            {
                temp[i*n+j] = A[row*n+col];
                j++;

                // Row is filled, so increase row index and
                // reset col index
                if (j == n - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
}

/* Recursive function for finding determinant of matrix.
   n is current dimension of A[][]. */
double determinant(double* A, int n,int m)
{
    double D = 0.0; // Initialize result

    //  Base case : if matrix contains single element
    if (m == 1)
        return A[0];

    //int temp[N][N]; // To store cofactors
    int i;
    double *temp = (double *)calloc(n*n, sizeof(double));

    int sign = 1;  // To store sign multiplier

     // Iterate for each element of first row
    int f;
    for (f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        getCofactor(A, temp, 0, f, n);
        D += sign * A[0*n+f] * determinant(temp,n, m - 1);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

// Function to get adjoint of A[N][N] in adj[N][N].
//void adjoint(int A[N][N],int adj[N][N])
void adjoint(double* A, double* adj, int n)
{
    if (n == 1)
    {
        adj[0] = 1;
        return;
    }

    // temp is used to store cofactors of A[][]
    int sign = 1;
    //temp[N][N];
    int i,j;
    double *temp = (double *)calloc(n*n, sizeof(double ));
    // for (i=0; i<n; i++)
    //      temp[i] = (double *)calloc(n, sizeof(double));


    for (i=0; i<n; i++)
    {
        for (j=0; j<n; j++)
        {
            // Get cofactor of A[i][j]
            getCofactor(A, temp, i, j, n);

            // sign of adj[j][i] positive if sum of row
            // and column indexes is even.
            sign = ((i+j)%2==0)? 1: -1;

            // Interchanging rows and columns to get the
            // transpose of the cofactor matrix
            adj[j*n+i] = (double)(sign)*(determinant(temp,n, n-1));
        }
    }

    free(temp);

}

// Function to calculate and store inverse, returns false if
// matrix is singular
void matInverse(double* A, double* inverse, int n)
{
    // Find determinant of A[][]
    double det = determinant(A,n,n);
    if (det == 0.0)
    {
        printf("Singular matrix, can't find its inverse");
        return;
    }
    // Find adjoint
    //int adj[N][N];
    int i,j;
    double *adj = (double *)calloc(n*n, sizeof(double));
    adjoint(A, adj, n);
    //display(adj,n);
    // Find Inverse using formula "inverse(A) = adj(A)/det(A)"
    for (i=0; i<n; i++)
        for (j=0; j<n; j++)
        {
            inverse[i*n+j] = (double)(adj[i*n+j]/det);
            printf("%f\n",inverse[i*n+j]);
        }

    free(adj);

    //return true;
}

void mydgetrf(double *arrA,int *pvt, int n, int b){
    int i,t,l,m,p,q,j,k,maxind,temps,end,ib;
    double tempv,max;
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
            //double *lln = (double*)calloc(sizeof(double), b*b);
            p=0;q=0;
            for(l=ib;l<=end;l++){
                for(m=ib;m<=end;m++){
                    if(l<m){
                        ll[p*b+q] = arrA[p*b+q];
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
            //matInverse(ll,lln,b);
            p=0;q=0;
            for(j=ib;j<=end;j++){
                //arrA[j*n+i] = arrA[j*n+i]/arrA[i*n+i];
                for(k=end+1;k<n;k++){
                    for(m=ib;m<end;m++){
                        arrA[j*n+k] += ll[p*b+q] * arrA[m*n+k];
                        q++;
                    }
                    q=0;
                }
                p++;
                q=0;
            }
        }
        for(j=end+1;j<n;j++){
            for(k=end+1;k<n;k++){
                for(l=ib;l<=end;l++){
                    arrA[j*n+k] = arrA[j*n+k] - arrA[j*n+l] * arrA[l*n+k];
                }
            }
        }
            free(ll);
            //free(lln);
        }
    }

}

void mydtrsm(int n, double *arrA, double *arrB, int *pvt, double *x, double *y, int label){
    double sum = 0.0, temp;
    int i,k;
    if(label == 0){
        y[0] = arrB[pvt[0]];
        for(i=1;i<n;i++){
            for(k=0;k<i-1;k++){
                sum+=y[k]*arrA[i*n+k];
            }
            y[i] = arrB[pvt[i]]-sum;
        }
    }
    else{
        x[n-1] = y[n-1]/arrA[n*n+n];
        for(i=n-2;i>=0;i--){
            for(k=i+1;k<n;k++){
                sum+= x[k]*arrA[i*n+k];
            }
            temp = y[i]-sum;
            x[i] = temp/arrA[i*n+i];

        }
    }
}

int main(){
    int *pvt,n,i,k;
    int ubound = 100, lbound = 0;
    int arrN[] = {1000,2000,3000,4000,5000};
    int block = 100;
    double *arrA, *arrB, *abk, *x, *y;
    double gflops, time;
    struct timespec tstart={0,0}, tend={0,0};
    int len = sizeof(arrN)/sizeof(arrN[0]);
    printf("\nBLOCKED GEPP \n");
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
        mydgetrf(arrA,pvt,n,block);
        clock_gettime(CLOCK_MONOTONIC, &tend);
        printf("\nSize = %d\n",n);

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
        printf("\n");
    }
    return 0;
}
