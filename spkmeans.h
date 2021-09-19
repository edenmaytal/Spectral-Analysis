
typedef struct cluster_element { /* implement each cluster as linked_list */
        double* vec; 
        struct cluster_element* next; 
    } cluster_element;

    typedef struct cluster_head { /* linked list of all clusters */
        struct cluster_element* centroid; 
        struct cluster_head* next; /* next cluster */
    } cluster_head;

typedef struct eigenData { 
        double eigenvalue; 
        double* eigenvector; 
        int order; /* for stable sorting */
} eigenData;

typedef struct vectors { 
        double** vecs; 
        int n; /* number of vectors */
} vectors;

typedef struct sp_result { 
        double** result; 
        int k; 
} sp_result; /* the final result, with the number of assigned clusters */

vectors* parseInput(char* file_name, int dim); /* return an array of vectors from input */
int get_dim_from_file(char* file_name);
sp_result* spectral_clustering(double** vecs, int dim, int n, int k, int is_py);
double l2_norm2vec(double* vec1, double* vec2, int dim); /* ||Xi-Xj|| */
double** wam(double** vecs, int dim, int n); /*Calculate weighted adjacency matrix */
double** ddg(double** wam, int n, int is_sqrt); /* caclculate ddg */
double** lnorm(double** ddg, double** wam,  int n); /* lnorm */
void calc_first_pivot(double** A, int n, int* pivot, int* max_index_array);
void update_pivot(double** A, int n, int* pivot,  int* max_index_array);
double** calcAprime(double** A, int n, int i, int j); 
double* calcIJ(double** A, int n); /* return a size 2 array where Aij is the max absoulte off-diagonal value */
double calc_theta(double** A, int i, int j);
double* eigen(double* A, int n); /* returns the diagonal values of A at the end of the iteration */
int find_k(eigenData** eigen, int n);
int stable_eigenvalues_sort(const void* a, const void* b);
double calc_t(double theta);
double calc_c(double t);
double calc_s(double t, double c);
double off2(double** A, int n);
void calc_Aprime(double** A, int n, int i, int j, double c, double s);
double** gen_identity_matrix(int n);
void update_V(double** V, int n, int i, int j, double c, double s);
double** gen_U(eigenData** V, int n, int k);
eigenData** jacobi(double** A, int n, int goal);
void gen_T( double** U, int n, int k);
int find_max_index(double** A, int n, int row);
int sign(double num);
void print_matrix(double** A, int n, int k);
double** kMeans(double** vecs,int K, int max_iter, int dim, int num_of_vecs);
double l2_norm(double* vec, double* centroid, int dimention);
void sum_vecs(double* sum_vec, double* vec, int d);
int eq_vecs(double* vec1, double* vec2, int d);
double* copyVec(double* vec, int dim);
cluster_head* init_clusters_list(double** vecs, int k);
cluster_head* find_closest_cluster(double* vec, cluster_head* head, int dim);
double* update_centroid(cluster_head* current_cluster ,int dim);
void empty_cluster(cluster_head* current_cluster, double* updated_centroid, int counter);
void print_matrix_columns(double** A, int n);
void print_vec(double* vec, int d);

