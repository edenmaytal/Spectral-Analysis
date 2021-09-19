#define PY_SSIZE_T_CLEAN
#include <Python.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "spkmeans.h"

static PyObject* fit(PyObject* self, PyObject* args);
static PyObject* fit_spk(PyObject* self, PyObject* args);
static PyObject* fit_kmeans(PyObject* self, PyObject* args);
static double* vecConvert(PyObject* vecPy, int dim);

static PyObject* fit(PyObject* self, PyObject* args)
{
    PyObject *_vecsPy, *vecPy;
    int dim, n, goal, i;
    Py_ssize_t n_py;
    double** vecs, **wam_mat, **ddg_res, **lnorm_mat;
    eigenData** eigen_array;
    if(!PyArg_ParseTuple(args, "iOi", &goal, &_vecsPy, &dim)){
        return NULL;
    }
    n_py= PyList_Size(_vecsPy); /*number of datapoints*/
    n = (int)n_py;
    vecs= (double**)calloc(n, sizeof(double*)); /*array of all datapoints*/
    if (vecs == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for(i=0; i<n; i++)
    {
        vecPy= PyList_GetItem(_vecsPy, i);
        vecs[i]=vecConvert(vecPy,dim); 
    }
    if (goal == 0) { /* wam */
        wam_mat = wam(vecs, dim, n);
        print_matrix(wam_mat, n, n);
        for (i=0; i<n; i++) {
            free(wam_mat[i]);
        }
        free(wam_mat); 
    }

    if (goal == 1) { /* ddg */
        wam_mat= wam(vecs, dim, n);
        ddg_res= ddg(wam_mat,n,0);
        print_matrix(ddg_res, n, n);
        for (i=0; i<n; i++) {
            free(wam_mat[i]);
        }
        free(wam_mat);
        for (i=0; i<n; i++) {
            free(ddg_res[i]);
        }
        free(ddg_res);
    }

    if (goal == 2) { /* lnorm */
        wam_mat= wam(vecs, dim, n);
        ddg_res= ddg(wam_mat,n,1);
        lnorm_mat= lnorm(ddg_res, wam_mat, n);
        print_matrix(lnorm_mat,n, n);
        for (i=0; i<n; i++) {
            free(lnorm_mat[i]);
        }
        free(lnorm_mat);
    }

    if (goal == 3) { /* jacobi */
        eigen_array = jacobi(vecs, n, 1);  
        free(eigen_array);      
    }
    Py_RETURN_NONE;
}

static PyObject* fit_spk(PyObject* self, PyObject* args)
{
    int k, dim, num_of_vecs, i, j;
    sp_result* result_struct;
    Py_ssize_t n;
    PyObject *_vecsPy, *lst, *vec, *num, *vecPy;
    double **vecs, **result;
    if(!PyArg_ParseTuple(args, "iOi", &k, &_vecsPy, &dim))
    {
        return NULL;
    }
    n= PyList_Size(_vecsPy); /*number of datapoints*/
    num_of_vecs = (int)n;
    vecs= (double**)calloc(num_of_vecs, sizeof(double*)); /*array of all datapoints*/
    if (vecs == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for(i=0; i<num_of_vecs; i++)
    {
        vecPy= PyList_GetItem(_vecsPy, i);
        vecs[i]=vecConvert(vecPy,dim); 
    }

    result_struct = spectral_clustering(vecs, dim, num_of_vecs, k, 1);
    dim = result_struct->k;
    result = result_struct->result;
    lst = PyList_New(n); /*convert the result array to a python list*/
    if (!lst || !PyList_Check(lst))
    {
        return NULL;
    }
    for (i = 0; i < num_of_vecs; i++)
    {
        vec = PyList_New(dim);
        if(!vec || !PyList_Check(vec))
        {
            return NULL;
        }
        for(j = 0; j < dim; j++)
        {
            num = PyFloat_FromDouble(result[i][j]); /*num= python(res[i][j]) */
            assert(type(num) == PyObject*);
            if(!num)
            {
                Py_DECREF(lst);
                return NULL;
            } 
            PyList_SetItem(vec, j, num); /*vec[j]=num */  
        }
        PyList_SetItem(lst, i, vec); /*lst[i]=vec*/
    }
    for (i=0; i<n; i++) {
        free(result_struct->result[i]);
    }
    free(result_struct->result);
    free(result_struct);
    return lst;
}

static PyObject* fit_kmeans(PyObject* self, PyObject* args) {
    PyObject *_vecsPy, *vecPy;
    int k, dim, num_of_vecs;
    Py_ssize_t i, n;
    double **vecs, **res;
    if(!PyArg_ParseTuple(args, "Oii", &_vecsPy, &k, &dim)) {
        return NULL;
    }
    n = PyList_Size(_vecsPy); /*number of datapoints*/
    vecs = (double**)calloc(n, sizeof(double*)); /*array of all datapoints*/
    if (vecs == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for(i=0; i<n; i++)
    {
        vecPy = PyList_GetItem(_vecsPy, i);
        vecs[i] = vecConvert(vecPy,dim); 
    }
    num_of_vecs = (int)n;
    res = kMeans(vecs, k, 300, dim, num_of_vecs); /*array of chosen centroids*/
    print_matrix(res ,k, k);
    for (i=0; i<k; i++) {
        free(res[i]);
    }
    free(res);
    Py_RETURN_NONE;
}
/* convert PyObject vec to C array*/
static double* vecConvert(PyObject* vecPy, int dim)
{
    Py_ssize_t i;
    PyObject* item;
    double* vec= (double*)calloc(dim, sizeof(double));
    if (vec == NULL) {
        printf("An Error Has Occured");
        exit(1);
    }
    for(i=0; i<dim; i++)
    {
        item= PyList_GetItem(vecPy, i);
        assert(type(item) == PyObject*);
        if (!PyFloat_Check(item)) {
            printf("An Error Has Occured");
        }
        vec[i]= PyFloat_AsDouble(item); 
        if (vec[i]  == -1 && PyErr_Occurred()){
            printf("An Error Has Occured"); 
        }        
    }
    return vec;
}

static PyMethodDef _methods[]= {
    {"fit", (PyCFunction) fit, METH_VARARGS, PyDoc_STR("Running all goals but 'spk'")},
    {"fit_spk", (PyCFunction) fit_spk, METH_VARARGS, PyDoc_STR("Running spk, returns matrix T")},
    {"fit_kmeans", (PyCFunction) fit_kmeans, METH_VARARGS, PyDoc_STR("Running kMeans")},
    {NULL, NULL, 0, NULL}
}; 

static struct PyModuleDef _moduledef = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    _methods
};

PyMODINIT_FUNC
PyInit_mykmeanssp(void)
{
    PyObject *m;
    m = PyModule_Create(&_moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}