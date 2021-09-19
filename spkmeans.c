#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "spkmeans.h"
#include <math.h>
#include <assert.h>
#define epsilon pow(10, -15)

/* get_dim_from_file: finds vector's dimention from file
input: 
file_name - name of file (extention .txt or .csv)
output: 
dim - datapoint's dimention */
int get_dim_from_file(char* file_name) {
    FILE* input;
    int i, dim = 0;
    char c;
    double coordinate;
    input = fopen(file_name, "r");
    if (input == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    i = fscanf(input, "%lf%c",&coordinate, &c);
        if (i == 0) {
            printf("An Error Has Occured");
            exit(1);
        }
        while (c != '\n'){
            if (c == ',') {
                dim++;
            }
            c = fgetc(input);
        }

    dim++;
    rewind(input);
    fclose(input);
    return dim;
} 
/* parseInput: parse input points from file to 2D array
input: 
file_name - name of file (extention .txt or .csv)
dim - datapoint's dimention 
output: 
datapoints - pointer to a struct (vectors*) contains: 
    vecs - array of datapoints (2D double array)
    n - number of datapoints  */
vectors* parseInput(char* file_name, int dim) {
    double coordinate;
    char c;
    int i=0, ind=0, d=0, is_first = 1;
    double* vec;
    FILE* input;
    vectors* datapoints;
    double** vecs = (double**)calloc(100, sizeof(double*)); 
     if (vecs == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    vec = (double*)calloc(dim, sizeof(double));
    if (vec == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    datapoints = (vectors*)malloc(sizeof(vectors));
     if (datapoints == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    input = fopen(file_name, "r");
     if (input == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }

    if (dim == 1) {
        while (fscanf(input, "%lf%c", &coordinate, &c) == 2) {
            if (is_first) {
                    is_first = 0;
                }
                else {
                    vec = (double*)calloc(dim, sizeof(double));
                    if (vec == NULL) {
                        printf("An Error Has Occured");
                        exit(1);
                    }
                }
            vec[0] = coordinate;
            vecs[ind] = vec;
            ind++;
        }
        vec = (double*)calloc(dim, sizeof(double));
            if (vec == NULL) {
                printf("An Error Has Occured");
                exit(1);
            }
        vec[i] = coordinate;
        vecs[ind]= vec; 
        ind++;
    }
    else {
    while (fscanf(input, "%lf%c", &coordinate, &c) == 2) { /* loop until EOF */
        d++;
        if (d != dim) {
            if (i == 0) {
                if (is_first) {
                    is_first = 0;
                }
                else {
                    vec = (double*)calloc(dim, sizeof(double));
                    if (vec == NULL) {
                        printf("An Error Has Occured");
                        exit(1);
                    }
                }
            }
            vec[i] = coordinate;
            i++;
        }
        else {
            d = 0;
            vec[i] = coordinate;
            i = 0;
            vecs[ind]= vec; 
            ind++;
        }
    }
    if (i != 0) {
        vec[i] = coordinate;
        vecs[ind]= vec; 
        ind++;
    }
    }
    vecs = realloc(vecs, ind*sizeof(double*));
    if (vecs == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    datapoints->vecs = vecs;
    datapoints->n = ind;
    fclose(input);
    return datapoints;
}

/* l2_norm2vec: perform l2 norm with 2 vectors 
inputs:
vec1, vec2 - double arrays 
dim - dimention of the vectors
output:
l2 norm of the two input vectors */
double l2_norm2vec(double* vec1, double* vec2, int dim)
{
    double norm;
    int i;
    norm=0; 
    for(i=0; i<dim; i++)
    {
        norm+=pow(vec1[i]-vec2[i],2);
    }
    return sqrt(norm);
}
/* wam: Calculationg wam for the given set of vectors
inputs:
vecs - 2D array of vectors
dim - the dimention of the vectors
n - number of vectors
output:
res - symmetric weighted adjacency matrix W nXn calculated from the vectors */
double** wam(double** vecs, int dim, int n)
{
    double** res;
    int i, j, curr_index=0;
    double norm;
    res= (double**)calloc(n, sizeof(double*));
    if (res == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for(i=0; i< n; i++)
    {
        res[i]= (double*)calloc(n, sizeof(double));
        if (res[i] == NULL) {
                printf("An Error Has Occured");
                exit(1);
            }
        for(j=0; j<n; j++)
            {
            if(i!=j)
            {
                norm= l2_norm2vec(*(vecs+i),*(vecs+j),dim);
                res[i][j]= exp((-1*norm/2));
                curr_index++;
            }
            else
            {
                res[i][j]=0;  /* no self loops */
            }
        }
    }
    for (i = 0; i < n; i++) {
        free(vecs[i]);
    }
    free(vecs);
    return res;
}
/* ddg: generate D or D^(-0.5) matrix from W matrix 
input: 
wam - symmetric weighted adjacency matrix W nXn (2D double array)
n - dimention of W matrix
is_squart - 0 if goal = 'ddg', else 1
output: 
res - if is_squart = 0 then res is diagonal degree matrix D nXn
      else res is D^(-0.5) */
double** ddg(double** wam, int n, int is_sqrt)
{
    int i, j; 
    double** res;
    double d=0; 
    res= (double**)calloc(n, sizeof(double*)); /* a n*n array with just 0s */
    if (res == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for(i=0; i<n; i++)
    {
        res[i]= (double*)calloc(n, sizeof(double));
        if (res[i] == NULL) {
            printf("An Error Has Occured");
            exit(1);
        }
        
        for(j=0; j<n; j++)
        {
            d+=wam[i][j];
        }
        if(is_sqrt==1)
        {
            res[i][i]= 1/sqrt(d);
        }
        else
        {
            res[i][i]= d;
        }
        d=0;
    } 
    return res;
}

/* lnorm: generate lnorm matrix from D^(-0.5) and W
input: 
ddg - D^(-0.5), D is diagonal degree matrix nXn
wam - symmetric weighted adjacency matrix W nXn (2D double array)
n - dimention of W and D matrixes
output: 
lnorm_matrix - normalized graph laplacian Lnorm nXn symmetric matrix */
double** lnorm(double** ddg, double** wam,  int n) {
    int i, j;
    double** lnorm_matrix = (double**)calloc(n, sizeof(double));
    if (lnorm_matrix == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for (i = 0 ; i < n; i++)
    {
        lnorm_matrix[i] = (double*)calloc(n, sizeof(double));
        if (lnorm_matrix[i] == NULL) {
            printf("An Error Has Occured");
            exit(1);
        }
    }
    for (i = 0; i < n; i++) {
        for (j = i; j < n; j++) {
            if (i == j) {
                lnorm_matrix[i][j] = 1;
            }
            else {
                lnorm_matrix[i][j] = (-1) * wam[i][j] * ddg[i][i] * ddg[j][j];
                lnorm_matrix[j][i] = lnorm_matrix[i][j];
            }
        }
    }
    for (i=0; i<n; i++) {
        free(wam[i]);
    }
    free(wam);
    for (i=0; i<n; i++) {
       free(ddg[i]);
    }
    free(ddg);
    return lnorm_matrix;
}
/* calc_first_pivot: calculates the first pivot of matrix A in Jacobi's algorythm
inout:
A -  symetric square matrix (if goal = 'spk' then A = Lnorm matrix)
n - dimention of A
pivot - int array
max_index_array - int array
output is void
pivot will contain the i,j indexes of matrix A's pivot
max_index_array - each cell i will contain the index of the max value in row i
 */
void calc_first_pivot(double** A, int n, int* pivot, int* max_index_array) { 
    int row, max_index_matrix;
    for (row = 0; row < n-1; row++) { /* go over A's rows*/
        max_index_array[row] = find_max_index(A, n, row);
    }
    /* find max index in the matrix*/
    max_index_matrix = 0;
    for (row = 0; row < n-1; row++) {
        if (fabs(A[row][max_index_array[row]]) > fabs(A[max_index_matrix][max_index_array[max_index_matrix]])) {
            max_index_matrix = row;
        }
    }
    pivot[0] = max_index_matrix;
    pivot[1] = max_index_array[max_index_matrix];
}

/*
update_pivot: calculates the pivot of matrix A in Jacobi's algorythm's
(not include the first iteration)
inputs: 
A - symetric square matrix (if goal = 'spk' then A = Lnorm matrix)
n - dimention of A matrix
pivot - int array [i,j], i,j are the indexes of the previous pivot  
max_index_array - each cell i will contain the index of the prev max value in row i
output is void
pivot will cointain pivot indexes i,j
max_index_array - each cell i will contain the index of the max value in row i
This function is *not* for the first iteration of Jacobi's algorythm 
*/
void update_pivot(double** A, int n, int* pivot,  int* max_index_array) {
    int row, i, j, max_index_matrix;
    i = pivot[0];
    j = pivot[1];
    for (row = 0; row < n-1; row++) { /* go over max index array */
        if ((max_index_array[row] = i) || (max_index_array[row] = j) || (row = 0) || (row = j)) { /* only i and j columns have changed */
            max_index_array[row] = find_max_index(A, n, row); /* find new row max index */
        }
        else {
            if (fabs(A[row][max_index_array[row]]) < fabs(A[row][i])) {
                max_index_array[row] = i;
            }
            if (fabs(A[row][max_index_array[row]]) < fabs(A[row][j])) {
                max_index_array[row] = j;
            }
        }
    }
     /* find max index in the matrix */
    max_index_matrix = 0;
    for (row = 0; row < n-1; row++) {
        if (fabs(A[row][max_index_array[row]]) > fabs(A[max_index_matrix][max_index_array[max_index_matrix]])) {
            max_index_matrix = row;
        }
    }
    pivot[0] = max_index_matrix;
    pivot[1] = max_index_array[max_index_matrix];
}
/*
find_max_index: finds the index of the row's abs max value in matrix A
inputs: 
A - symetric square matrix
n - dimention of A matrix
row - current row index
output: input of the row's max abs value (off diagonal)
 */
int find_max_index(double** A, int n, int row) {
    int col;
    int max_index_row = row + 1; /* first off-diag index */
        for (col = row + 1; col < n; col++) { /*go over elements off-diagonal (upper half) */
            if (fabs(A[row][col]) > fabs(A[row][max_index_row])) { 
                max_index_row = col;
            }
        }
    return max_index_row;
}
/*
input: 
A - symetric square matrix
i,j - pivot indices
output: theta  
*/
double calc_theta(double** A, int i, int j) {
    double theta;
    theta = (A[j][j] - A[i][i]) / (2*A[i][j]);
    return theta;
}
/*
input: theta 
output: t  
*/
double calc_t(double theta) {
    int s = sign(theta);
    double t;
    t = s / (fabs(theta) + sqrt(pow(theta,2)+1));
    return t;
}
/*
input: num - a double 
output: 1 if the number is positive or 0, else -1 
*/
int sign(double num) {
    if (num >= 0) return 1;
    else return -1;
}
/*
input: t 
output: c  
*/
double calc_c(double t) {
    double c = 1 / sqrt(pow(t,2) + 1);
    return c; 
}
/*
input: t, c 
output: s  
*/
double calc_s(double t, double c) {
    double s = t*c;
    return s;
}
/* calc_Aprime: calc A prime (symetric square matrix) in place
inputs: 
A - symetric square matrix
n - dimention of A matrix 
i,j - pivot indices
c,s
output is void
 */
void calc_Aprime(double** A, int n, int i, int j, double c, double s) {
    int r, ind;
    double* I; /*values of A[][i]*/
    double* J; /*values of A[][j]*/
    I = (double*)calloc(n, sizeof(double));
    if (I == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    J = (double*)calloc(n, sizeof(double));
    if (J == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for (ind = 0; ind < n; ind++) {
        I[ind] = A[ind][i];
    }
     for (ind = 0; ind < n; ind++) {
        J[ind] = A[ind][j];
    }
    for (r = 0; r < n; r++) {
        if ((r != i) && (r != j)) {
            A[r][i] = c*I[r] - s*J[r]; /*c*A[r][i] - s*A[r][j] */
            A[i][r] = c*I[r] - s*J[r];
            A[r][j] = c*J[r] + s*I[r]; /* c*A[r][j] + s*A[r][i] */
            A[j][r] = c*J[r] + s*I[r];
        } 
    }
    A[i][i] = (pow(c,2)*I[i]) + (pow(s,2)*J[j]) - (2*s*c*J[i]);
    A[j][j] = (pow(s,2)*I[i]) + (pow(c,2)*J[j]) + (2*s*c*J[i]);
    A[i][j] = 0;
    A[j][i] = 0;
    free(I);
    free(J);
}

/* off2: calculte off^2(A) - sum of squares of all off-diagonal elements of A
input:
A - symmetric square matrix
n - dimention of A
output:
sum_sq_off_diag - off^2 value of A matrix 
*/
double off2(double** A, int n) {
    double sum_sq_off_diag = 0;
    int i,j;
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++) {
            if (i != j) { 
                sum_sq_off_diag += pow(A[i][j], 2);
            }
        }
    }
    return sum_sq_off_diag;
} 

/* gen_identity_matrix: generate identinty matrix 
input: n - dimention of the matrix
output: identity nXn matrix
*/
double** gen_identity_matrix(int n) {
    double *p;
    double **I;
    int i;
    p = (double*)calloc(n*n, sizeof(double));
    if (p == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    I = (double**)calloc(n, sizeof(double*));
    if (I == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for (i = 0; i < n; i++) {
        I[i] = p + (i*n);
    }
    for (i = 0; i < n; i++) {
        I[i][i] = 1;
    }
    return I;
}
/* update_V: update V matrix (eigen vectors matrix) in place
inputs: 
V - previous V matrix
n - dimention of V matrix 
i,j - pivot indices
c,s
output is void
 */
void update_V(double** V, int n, int i, int j, double c, double s) {
    int r, ind;
    double* I; /*values of V[][i]*/
    double* J; /*values of V[][j]*/
    I = (double*)calloc(n, sizeof(double));
    if (I == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    J = (double*)calloc(n, sizeof(double));
    if (J == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for (ind = 0; ind < n; ind++) {
        I[ind] = V[ind][i];
    }
     for (ind = 0; ind < n; ind++) {
        J[ind] = V[ind][j];
    }
    for (r = 0; r < n; r++) {
        V[r][i] = c*I[r] - s*J[r]; /* c*V[r][i] - s*V[r][j] */
        V[r][j] = s*I[r] + c*J[r]; /* s*V[r][i] + c*V[r][j] */
    }
    free(I);
    free(J);
}
/*
jacobi: perform jacobi's algorythm
calculates eigenvalues and eigenvectors of symmetric matrix A
input:
A - symetric square matrix
n - dimention of the matrix
goal - 1 if the goal is 'jacobi', else 0
outupt:
eigen_array - array of pointer to eigenData struct
eigenData contain eigenvalue, matching eigenvale and array index 
*/
eigenData** jacobi(double** A, int n, int goal) {
    int k, iter = 1;
    int i,j;
    double** V = gen_identity_matrix(n);
    double* vec;
    double theta, t,c,s, offA2, offA_prime;
    int pivot_indexes[2]; 
    eigenData** eigen_array;
    eigenData* eigen;
    int* max_index_array = (int*)calloc(n-1, sizeof(int));
    if (max_index_array == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    if (max_index_array == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    eigen_array = (eigenData**)calloc(n, sizeof(eigenData*));
    if (eigen_array == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    calc_first_pivot(A, n, pivot_indexes, max_index_array);
    theta = calc_theta(A, pivot_indexes[0], pivot_indexes[1]);
    t = calc_t(theta);
    c = calc_c(t);
    s = calc_s(t,c);
    offA2 = off2(A, n);
    i = pivot_indexes[0];
    j = pivot_indexes[1];
    calc_Aprime(A, n, i, j, c, s);
    offA_prime = off2(A, n);
    update_V(V, n, i, j, c, s);
        
    while (((offA2 - offA_prime) > epsilon) && (iter < 100)) {
        update_pivot(A, n, pivot_indexes, max_index_array);
        theta = calc_theta(A, pivot_indexes[0], pivot_indexes[1]);
        t = calc_t(theta);
        c = calc_c(t);
        s = calc_s(t,c);
        i = pivot_indexes[0];
        j = pivot_indexes[1];
        update_V(V, n, i, j, c, s);
        offA2 = offA_prime;
        calc_Aprime(A, n, i, j, c, s);
        offA_prime = off2(A, n);
        iter++;
    }
    if (goal) {
        for (k = 0; k < n-1; k++) { /*print eigenvalues */
            if ((A[k][k] < 0) && (A[k][k] > -0.00005)) {
                A[k][k] = 0;
            }
            printf("%.4f%s", A[k][k], ",");
        }
        if ((A[n-1][n-1] < 0) && (A[n-1][n-1] > -0.00005)) {
            A[n-1][n-1] = 0;
        }
        printf("%.4f%s", A[n-1][n-1], "\n");
        print_matrix_columns(V,n); /*print eigenvectors */
    }
    
    else {
        for (k = 0; k < n; k++) {
            eigen = (eigenData*)malloc(sizeof(eigenData));
            if (eigen == NULL) {
                printf("An Error Has Occured");
                exit(1);
            }
            eigen->eigenvalue= A[k][k];
            eigen->order = k;
            vec=(double*)calloc(n, sizeof(double));
            if (vec == NULL) {
                printf("An Error Has Occured");
                exit(1);
            }
            for(i = 0; i < n; i++) 
            {
            vec[i]= V[i][k];
            }
            eigen->eigenvector= vec;
            eigen_array[k] = eigen;
        }
    }
    free(max_index_array);
    free(V[0]); 
    free(V);
    for (i = 0; i < n; i++) {
        free(A[i]);
    }
    free(A);
    return eigen_array;
}
/*
find_k: find k using eigengap Heuristic
input:
eigen - array of pointer to eigenData struct 
        each struct contain eigenvalue, eigenvector and index
n - number of pointers in eigen
ouptup: k
*/
int find_k(eigenData** eigen, int n) {
    int i, max_gap_index = 0;
    double max_gap = 0, current_gap = 0;
    qsort(eigen, n, sizeof(eigenData*), stable_eigenvalues_sort);
    for (i = 0; i < n/2; i++) { 
        current_gap = fabs(eigen[i]->eigenvalue - eigen[i+1]->eigenvalue); /*delta_i */
        if (current_gap > max_gap) {
            max_gap = current_gap;
            max_gap_index = i;
        }
    }
    return  max_gap_index + 1;
}
/* inputs:
a,b - pointers to eigenData struct elements
output:
-1 if a's eigenvalue < b's eigenvalue, else 1
if a's eigenvalue = b's eigenvalue, the result will be according to their index (order)
*/
int stable_eigenvalues_sort(const void* a, const void* b) {
    const eigenData* e1 = *(const eigenData **)a;
    const eigenData* e2 = *(const eigenData **)b;
    if (e1->eigenvalue < e2->eigenvalue) {
        return -1;
    }
    else if (e1->eigenvalue > e2->eigenvalue) {
        return 1;
    }
    else if (e1->order < e2->order) {
        return -1;
    }
    else if (e1->order > e2->order) {
        return 1;
    }
    else {
        return 0;
    }
}

/* print_matrix: prints matrix nXk
input:
A - matrix to print (2D double array)
n - first dimention of A
k - second dimention of A
 */
void print_matrix(double** A, int n, int k) {
    int row, col;
    for (row = 0; row < n-1; row++) {
        for (col = 0; col < k-1; col++) {
            if ((A[row][col] < 0) && (A[row][col] > -0.00005)) {
                A[row][col] = 0;
            }
            printf("%.4f%s", A[row][col], ",");
        }
        if ((A[row][k-1] < 0) && (A[row][k-1] > -0.00005)) {
            A[row][k-1] = 0;
        }
        printf("%.4f%s", A[row][k-1], " \n");
    }
    for (col = 0; col < k-1; col++) {
        if ((A[row][col] < 0) && (A[row][col] > -0.00005)) {
            A[row][col] = 0;
        }
        printf("%.4f%s", A[row][col], ",");
    }
     if ((A[row][k-1] < 0) && (A[row][k-1] > -0.00005)) {
            A[row][k-1] = 0;
        }
        printf("%.4f", A[row][k-1]);
}
/* print_matrix_columns: prints square matrix nXn, each printed row is a column in the matrix
ipnut:
A - matrix to print (2D double array)
n - dimention of A
 */
void print_matrix_columns(double** A, int n) {
    int row, col;
    for (col = 0; col < n-1; col++) {
        for (row = 0; row < n-1; row++) {
            if ((A[row][col] < 0) && (A[row][col] > -0.00005)) {
                A[row][col] = 0;
            }
            printf("%.4f%s", A[row][col], ",");
        }
        if ((A[n-1][col] < 0) && (A[n-1][col] > -0.00005)) {
            A[n-1][col] = 0;
        }
        printf("%.4f%s", A[n-1][col], " \n");
    }
    for (row = 0; row < n-1; row++) {
            if ((A[row][col] < 0) && (A[row][col] > -0.00005)) {
                A[row][col] = 0;
            }
            printf("%.4f%s", A[row][col], ",");
        }
        if ((A[n-1][col] < 0) && (A[n-1][col] > -0.00005)) {
            A[n-1][col] = 0;
        }
        printf("%.4f", A[n-1][col]);
}
/*generate U: nXk matrix with eigenvectors as columns
inputs:
V - array of pointers to eigenData struct (ordered by the eigenvalues)
n - number of eigenvecrots in V
k - number of vectors as columns of U
output:
U - matrix U
*/
double** gen_U(eigenData** V, int n, int k) 
{
    double** U;
    double* matrix_row;
    int row, col, i; 
    U = (double**)calloc(n, sizeof(double*));
    if (U == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for (row = 0; row < n; row++) {
        matrix_row = (double*)calloc(k, sizeof(double));
        if (matrix_row == NULL) {
            printf("An Error Has Occured");
            exit(1);
        }
        if (matrix_row == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
            for(col = 0; col < k; col++) 
            {
                matrix_row[col]= (V[col])->eigenvector[row];
            }
        U[row] = matrix_row;
    }
    for (i = 0; i<n; i++) {
        free(V[i]->eigenvector);
        free(V[i]);
    }
    free(V);
    return U;
}
/* gen_T: generate matrix T by normalize every row in matrix U (in place)
inputs:
U - matrix nXk
n - first dimention of U
k - second dimention of U
output is void, matrix T is created in place from matrix U */
void gen_T( double** U, int n, int k)
{
    int i,j;
    double norm;
    double* zero_vec = (double*)calloc(n, sizeof(double));
    if (zero_vec == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for(i = 0; i < n; i++)
    {
        norm = l2_norm2vec(U[i], zero_vec, k);
        if (norm == 0) {
            norm = 1;
        }
        for(j=0; j<k; j++)
        {
            U[i][j]= U[i][j]/norm;
        }
    }
    free(zero_vec); 
}

/* spectral_clustering: main algorythm to perform spectral clustering
inputs:
vecs - the vectors to cluster (2D double array)
dim - vector's dimension 
n - number of the vectors
k - number of clusters (0 if not specified)
is_py - 1 if the function was called from python, else 0
output:
result_struct - pointer to a struct (sp_result) contains:
    T - array of double array
    k - updated number of clusters (if the input was k = 0)
if is_py = 0, prints the final centroids after perform kMean */

sp_result* spectral_clustering(double** vecs, int dim, int n, int k, int is_py)
{
    int tmp, i = 0;
    double** wam_mat, **ddg_res, **lnorm_mat, **T, **result;
    eigenData** eigen_array;
    sp_result* result_struct = (sp_result*)malloc(sizeof(sp_result));
    if (result_struct == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    wam_mat= wam(vecs, dim, n);
    ddg_res= ddg(wam_mat,n,1);
    lnorm_mat= lnorm(ddg_res, wam_mat, n);
    eigen_array = jacobi(lnorm_mat, n, 0);
    tmp = find_k(eigen_array,n);
    if(k==0)
    {
        k= tmp;
    }
    T=gen_U(eigen_array,n, k);
    gen_T(T,n,k);
    result_struct->result = T;
    result_struct->k = k;
    if (!is_py) { /* result will be the final centroids */
        result = kMeans(T, k, 300, k, n);
        print_matrix(result ,k, k);
        for (i = 0; i < k; i++) {
            free(result[i]);
        }
        free(result);
    }
    return result_struct;
}

/* kMeans Algorythm */

/* init_clusters_list: initiate linked list of clusters from datapoints 2D array
inputs:
vecs - 2D double array of datapoints
k - number of clusters in the list
output:
head - pointer to the head of the list (cluster_head struct)
 */ 
cluster_head* init_clusters_list(double** vecs, int k) {
    int i;
    cluster_head *head, *new_cluster;
    cluster_element *current_element;
    cluster_head *current_cluster = (cluster_head*)malloc(sizeof(cluster_head));
    if (current_cluster == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    head = current_cluster; /* pointer to the beginning of the clusters list */
    for (i = 0; i < k-1; i++) { /* initialize all clusters */
        current_element = (cluster_element*)malloc(sizeof(cluster_element));
        if (current_element == NULL) {
            printf("An Error Has Occured");
            exit(1);
        }
        current_element ->vec = *(vecs+i);
        current_element-> next = NULL;
        current_cluster->centroid = current_element;

        new_cluster = (cluster_head*)malloc(sizeof(cluster_head));
        if (new_cluster == NULL) {
            printf("An Error Has Occured");
            exit(1);
        }
        current_cluster->next = new_cluster;
        current_cluster = new_cluster;
    }
    current_element = (cluster_element*)malloc(sizeof(cluster_element));
    if (current_element == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    current_element ->vec = *(vecs+i);
    current_element-> next = NULL;
    current_cluster->centroid = current_element;
    current_cluster->next = NULL;
    return head;
}

/* find_closest_cluster: 
input: 
vec - vector (double array) for associating to the closest cluster
head - pointer to the head of the cluster's linked list
dim - dimention of vec
output:
p_min - pointer to the closest cluster (cluster_head struct) */    
cluster_head* find_closest_cluster(double* vec, cluster_head* head, int dim) {
    cluster_head *current_cluster = head;
    cluster_element* current_centroid = current_cluster-> centroid; 
    double min_dist = l2_norm(vec, current_centroid-> vec, dim); 
    double dist;
    cluster_head *p_min = head; /*pointer to the current closest cluster  */
    current_cluster = current_cluster ->next;
    while (current_cluster != NULL) { /* go over all cluster */
            current_centroid = current_cluster -> centroid;
            dist = l2_norm(vec, current_centroid -> vec, dim);
            if (dist < min_dist) { 
                p_min = current_cluster;
                min_dist = dist;
            }
            current_cluster = current_cluster-> next;               
        }
    return p_min;
}

/* update_centroid:
input: 
current_cluster - pointer to cluster (cluster_head struct) 
dim - dimention of cluster's vector
output:
sum_vec - value of current_cluster's updated centroid*/
double* update_centroid(cluster_head* current_cluster ,int dim) {
    double *sum_vec;
    int i, cluster_len = 0;
    cluster_element* curr = current_cluster->centroid->next; /*first element in cluster*/
    sum_vec = (double*)calloc(dim, sizeof(double));
    if (sum_vec == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }  
    while(curr!=NULL) { /*go over cluster vevtors */
        sum_vecs(sum_vec, curr->vec, dim); 
        cluster_len+=1;          
        curr=curr->next;    
    }          
    for(i = 0; i < dim; i++) {
        if(cluster_len!=0) {
        sum_vec[i] = sum_vec[i]/cluster_len;
        }
    }
    return sum_vec;
}

/* empty_cluster: remove old vectors from cluster (after updating the centroid)
input:
current_cluster - pointer to cluster (cluster_head struct)
updated_centroid - new centroid value 
counter 
output is void
 */
void empty_cluster(cluster_head* current_cluster, double* updated_centroid, int counter) {
    double* vec_to_free = current_cluster->centroid->vec;
    cluster_element* element_to_free, *tmp;
    current_cluster->centroid->vec = updated_centroid;
    if (counter != 1) {
        free(vec_to_free);
    } 
    element_to_free = current_cluster->centroid->next;
    current_cluster-> centroid->next= NULL;
    if (element_to_free != NULL) {
        while (element_to_free != NULL) {
            tmp = element_to_free;
            element_to_free = element_to_free-> next;
            free(tmp);
        }
    } 
}
/* kMeans: perform clustring using kMeans algorythm
inputs:
vecs - 2D double array of datapoints to cluster
k - required number of clusters
max_iter - max number of iterations for the algorythm
dim - dimention of datapoints
num_of_vecs - number to datapoints in vecs
output:
res - 2D double array of k final centroids
 */      
double** kMeans(double** vecs, int k, int max_iter, int dim, int num_of_vecs) {
    double** res;
    int i, counter = 0; 
    int is_convergence = 0, is_equal_centroid; /*booleans*/
    double *updated_centroid;
    cluster_head *current_cluster, *head, *p_min, *tmp1;
    cluster_element *current_centroid, *new_element;
    head = init_clusters_list(vecs, k);
    while ((counter < max_iter + 1) && !is_convergence) {
        counter++;
        for (i = 0; i < num_of_vecs; i++) { /* go over all vectors xi */
            p_min = find_closest_cluster(vecs[i], head, dim);
            /* assign vec xi to the closest cluster p_min */
            new_element = (cluster_element*)malloc(sizeof(cluster_element)); 
            if (new_element == NULL) {
                printf("An Error Has Occured");
                exit(1);
            }
            new_element -> vec = vecs[i];
            current_centroid = p_min-> centroid; /* the first element of the closest cluster */
            if (current_centroid-> next == NULL) {
                current_centroid-> next = new_element;
                new_element-> next = NULL;
            } 
            else {
                new_element -> next = current_centroid-> next;
                current_centroid-> next = new_element;
            }
        }
        /* update centroid */
        is_convergence = 1;
        current_cluster = head;
        while (current_cluster!= NULL) { /* go over clusters list*/
            updated_centroid = update_centroid(current_cluster, dim);
            is_equal_centroid = eq_vecs(updated_centroid, current_cluster->centroid->vec,dim);
            if (is_equal_centroid) {
                is_convergence = 0;
            }
            empty_cluster(current_cluster, updated_centroid, counter);
            current_cluster= current_cluster->next;
        }
    }
    current_cluster= head; 
    res= (double**)calloc(k, sizeof(double*)); /*array of final centroids to return*/
    if (res == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    i=0;
    while(current_cluster!=NULL) {
        current_centroid= current_cluster->centroid;
        res[i]=copyVec(current_centroid->vec,dim);
        free(current_centroid->vec);
        free(current_centroid);
        tmp1= current_cluster;
        current_cluster= current_cluster->next;
        i++; 
        free(tmp1);
    }
    for(i=0; i < num_of_vecs; i++) {
        if(*(vecs+i)!=NULL) {
            free(*(vecs+i));
        }
    }
    free(vecs);
    return res;
}

/* copyVec: 
inputs:
vec - vector to copy (double array)
dim - dimention of the vector
ouput:
copy of the vector (double array) */
 double* copyVec(double* vec, int dim)
{
    int i=0;
    double* res= (double*)calloc(dim, sizeof(double));
     if (res == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for(i=0; i<dim; i++)
    {
        res[i]= vec[i];
    }
    return res;
}
/* eq_vecs: return 1 if the two vectors are not equal, or 0 otherwise
inputs:
vec1, vec2 - vectors to comapre (double array)
dim - dimention of vec1 and vec2
output:
1 if vec1 = vec2, else 0 */ 
int eq_vecs(double* vec1, double* vec2, int dim)
{
   int i;
    for(i=0; i<dim; i++)
    {
        if(*(vec1+i)!=*(vec2+i))
        {
            return 1;
        }
    }
    return 0;
}
/* sum_vecs: sum two vectors element by element
input:
sum_vec - first vector to sum (will be the result vector)
vec - second vector to sum
dim - dimention of the vectors
output is void */
void sum_vecs(double* sum_vec, double* vec, int dim)
{
    int i;
    for(i=0; i<dim; i++)
    {
        *(sum_vec+i)+= *(vec+i);
    }
}
 /*l2_norm: calculate l2 norm of two vectors
 input:
 vec - 2D double array
 centroid - 2D double array
 output:
 sum1 - l2 norm value */
double l2_norm(double* vec, double* centroid, int dimention) {
    double sum1 = 0;
    int i;
    for (i=0; i < dimention; i++){
        sum1 += (vec[i] - centroid[i])*(vec[i] - centroid[i]);
    }
    return sum1;
}

int main(int argc, char** argv){
    int i;
    sp_result* result_struct;
    int k = atoi(argv[1]);
    int dim = get_dim_from_file(argv[3]);
    vectors* datapoints = parseInput(argv[3], dim);
    double** vecs = datapoints->vecs;
    int n = datapoints->n;
    double **wam_mat, **ddg_res, **lnorm_mat;
    eigenData** eigen_array;
    const char* goal = argv[2];
    if (argc != 4) {
        i = 0;
    }
    
    if (!strcmp(goal, "wam")) {
        wam_mat = wam(vecs, dim, n);
        print_matrix(wam_mat, n, n); 
        for (i=0; i<n; i++) {
            free(wam_mat[i]);
        }
        free(wam_mat);
        free(datapoints);
    }

    if (!strcmp(goal, "ddg")) {
        wam_mat= wam(vecs, dim, n);
        free(datapoints);
        ddg_res= ddg(wam_mat,n,0);
        for (i=0; i<n; i++) {
            free(wam_mat[i]);
        }
        free(wam_mat);
        print_matrix(ddg_res, n, n);
        for (i=0; i<n; i++) {
            free(ddg_res[i]);
        }
        free(ddg_res);
    }
    if (!strcmp(goal, "lnorm")) {
        wam_mat= wam(vecs, dim, n);
        free(datapoints);
        ddg_res= ddg(wam_mat,n,1);
        lnorm_mat= lnorm(ddg_res, wam_mat, n);
        print_matrix(lnorm_mat, n, n);
        for (i=0; i<n; i++) {
            free(lnorm_mat[i]);
        }
        free(lnorm_mat);    
    }

    if (!strcmp(goal, "jacobi")) {
        eigen_array = jacobi(vecs, n, 1);
        free(datapoints);
        free(eigen_array);
    }
    if (!strcmp(goal, "spk")) {
        result_struct = spectral_clustering(vecs, dim, n, k, 0);
        free(result_struct);
        free(datapoints);
    }
    return 0;
}

