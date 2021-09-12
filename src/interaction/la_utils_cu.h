/*! PArticle Wake ANalysis
 * \file la_utils_cu.h
 * \brief Inline functions for matrix data conversion from gsl library to regular C++ pointers
 * (gsl is not CUDA compatible)
 *
 * @author 	Sumeet Kumar
 * @date	09/07/2021
 *
 */

#ifndef LA_UTILS_CU_H
#define LA_UTILS_CU_H

#include <float.h>

#define SOFTENING       DBL_EPSILON

inline void la_alloc( double **uu, const size_t &row_size, const size_t &col_size){
    uu = new double* [row_size];
    for ( size_t row = 0; row < row_size; ++row ) {
        uu[row] = new double[col_size];
    }
};

inline double** la_calloc( const size_t &row_size, const size_t &col_size){
    double** uu = new double* [row_size];
    for ( size_t row = 0; row < row_size; ++row ) {
        uu[row] = new double[col_size];
        for (size_t col = 0; col < col_size; ++col) {
            uu[row][col] = 0.0;
        }
        return uu;
    }
};

/*!
 *allocates a matrix array and copies elements of gsl matrix to it
 * (every la_gslalloc() should be followed by la_gsldealloc())
 */
inline double** la_gslalloc( const gsl_matrix *vv, const size_t &col_size,  const size_t &row_size){
    double** uu =  new double* [col_size];
    for ( size_t row = 0; row < col_size; ++row ) {
        uu[row] = new double[row_size];
    }
    for ( size_t row = 0; row < col_size; ++row ) {
        for (size_t col = 0; col < row_size; ++col) {
            uu[row][col] = gsl_matrix_get(vv, row, col);
        }
    }
    return uu;
};

inline double* la_gslalloc( const gsl_vector *v,  const size_t &row_size){
    double* u =  new double [row_size];
    for (size_t col = 0; col < row_size; ++col) {
        u[col] = gsl_vector_get(v, col);
    }
    return u;
};

inline void la_gslalloc_reverse( double **uu, gsl_matrix *vv,  const size_t &row_size, const size_t &col_size){
    for ( size_t row = 0; row < row_size; ++row ) {
        for (size_t col = 0; col < col_size; ++col) {
            gsl_matrix_set(vv, row, col, uu[row][col]);
        }
    }
};
inline void la_gslalloc_reverse(const double *u, gsl_vector *v,  const size_t &row_size){
    for (size_t col = 0; col < row_size; ++col) {
        gsl_vector_set(v, col, u[col]);
    }
};

inline void la_dealloc( double **uu,  size_t &col_size){
    for (size_t row = 0; row < col_size; ++row) {
        delete[] uu[row];
    }
    delete[] uu;
};

inline void la_dealloc( double *u){
    delete[] u;
};
inline void la_gsl_cross(const double *a, const double *b, double *c, const size_t &row_size){//note that this works only for 3 dim vectors??
    if(row_size==3){
        double ax = a[0];
        double ay = a[1];
        double az = a[2];
        double bx = b[0];
        double by = b[1];
        double bz = b[2];
        c[0] = ay*bz-az*by;
        c[1] = az*bx-ax*bz;
        c[2] = ax*by-ay*bx;
    }
};

inline void la_gsl_blas_ddot(const double *u, const double *v, double *result, const size_t &row_size){
    for(size_t col = 0; col<row_size; ++col){
        *result = *result + u[col]*v[col];
    }
};

inline double la_gsl_blas_dnrm2_soft(const double *u,const size_t &row_size){
    double result = 0;
    for(size_t col = 0; col<row_size; ++col){
        result = result + u[col]*u[col];
    }
    result = sqrt(result+DBL_EPSILON);
    return result;
};

inline double la_gsl_blas_dnrm2(const double *u,const size_t &row_size){
    double result = 0;
    for(size_t col = 0; col<row_size; ++col){
        result = result + u[col]*u[col];
    }
    result = sqrt(result);
    return result;
};

inline void la_gsl_blas_dscal(const double &alpha, double *u, const size_t &row_size){
    for(size_t col = 0; col<row_size; ++col){
        u[col] = alpha*u[col];
    }
};

inline void la_gsl_vector_add( double *u,const double *v, const size_t &row_size){
    for(size_t col = 0; col<row_size; ++col){
        u[col] = u[col]+v[col];
    }
};

inline void la_gsl_vector_sub( double *u,const double *v, const size_t &row_size){
    for(size_t col = 0; col<row_size; ++col){
        u[col] = u[col]-v[col];
    }
};

inline double la_gsl_pow_2(const double w){
    return w*w;
};

inline void la_gsl_vector_memcpy( double *u,const double *v, const size_t &row_size){
    for(size_t col = 0; col<row_size; ++col){
        u[col] = v[col];
    }
};
/*
inline double* la_gsl_matrix_row(const double **uu,const size_t &row_num,const size_t &row_size){
    double *v = uu;
    for(size_t col = 0; col<row_size; ++col){
        v[col] = uu[row_num][col];
    }
    return v;
};
*/

inline double* la_gsl_matrix_row(double **uu,const size_t &row_num) {
    double *v = uu[row_num];
    return v;
}

inline double la_gsl_vector_get(const double *u,const size_t &col_num){
    return u[col_num];
};

inline void la_gsl_vector_set(double *u, const size_t &col_num, const double &v){
    u[col_num] = v;
};
/*
inline void la_gsl_vector_memcpy(double *u, double *v, const size_t &row_num){
    return u[row_num];
};
*/
inline void set_gsl_drow(double *u, const double **vv, size_t &row, size_t &col_size){
    for(size_t col = 0; col<col_size; ++col){
        u[col] = vv[row][col];
    }
};

inline void la_printVec(const double *uu, size_t &row_size){
    printf("using printVec");
    for(size_t col = 0; col<row_size; ++col) {
        printf("\t%f", uu[col]);
    }
    printf("\n");
};


inline void la_print( double **uu,  size_t &row_size, size_t &col_size){
    for ( size_t row = 0; row < row_size; ++row ) {
        for (size_t col = 0; col < col_size; ++col) {
            printf("\t%f", uu[row][col]);
        }
        printf("\n");
    }
    printf("\n");
};

inline void la_print(const double *u,  size_t &row_size){
    for (size_t col = 0; col < row_size; ++col) {
        printf("\t%.16f", u[col]);
    }
    printf("\n");
};


#endif
