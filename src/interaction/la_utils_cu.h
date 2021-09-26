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

#include <stdio.h>
#include <float.h>

#define SOFTENING       DBL_EPSILON

__host__ __device__
double4* la_to_d4alloc(double **uu, const size_t &col_size, const size_t &row_size){
    int bytes = col_size*sizeof(double4);
    double *buf = (double*)malloc(bytes);
    double4 *v = (double4*)buf;
    for ( size_t row = 0; row < col_size; ++row ) {
        if(row_size == 3) {
            v[row].x = uu[row][0];
            v[row].y = uu[row][1];
            v[row].z = uu[row][2];
            v[row].w = 0.0;
        }else if(row_size == 2){
            v[row].x = uu[row][0];
            v[row].y = uu[row][1];
            v[row].z = 0.0;
            v[row].w = 0.0;
        }
    }
    return v;
};

__host__ __device__
double4* la_to_d4alloc(double *u, const size_t &col_size){
    int bytes = col_size*sizeof(double4);
    double *buf = (double*)malloc(bytes);
    double4 *v = (double4*)buf;
    for ( size_t row = 0; row < col_size; ++row ) {
        v[row].x = u[row];
        v[row].y = 0.0;
        v[row].z = 0.0;
        v[row].w = 0.0;
    }
    return v;
};

__host__ __device__
double4* la_d4calloc(const size_t &col_size){
    int bytes = col_size*sizeof(double4);
    double *buf = (double*)malloc(bytes);
    double4 *v = (double4*)buf;
    for ( size_t row = 0; row < col_size; ++row ) {
            v[row].x = 0.0;
            v[row].y = 0.0;
            v[row].z = 0.0;
            v[row].w = 0.0;
    }
    return v;
};

double4* la_to_d4alloc_reverse(double **uu,double4 *v, const size_t &col_size, const size_t &row_size){
    for ( size_t row = 0; row < col_size; ++row ) {
        if(row_size == 3) {
            uu[row][0] = v[row].x;
            uu[row][1] = v[row].y;
            uu[row][2] = v[row].z;
        }else if(row_size == 2){
            uu[row][0] = v[row].x;
            uu[row][1] = v[row].y;
        }
    }
    return v;
};

__device__
void la_d4_cross(const double4 *a, const double4 *b, double4 *c){
    double ax = a->x;
    double ay = a->y;
    double az = a->z;
    double bx = b->x;
    double by = b->y;
    double bz = b->z;
    c->x = ay*bz-az*by;
    c->y = az*bx-ax*bz;
    c->z = ax*by-ay*bx;
    c->w = 0.0;
};

__device__
void la_d4_blas_ddot(const double4 *u, const double4 *v, double *result){
    *result = u->x*v->x + u->y*v->y + u->z*v->z + u->w*v->w;
};

__device__
double la_d4_blas_dnrm2_soft(const double4 *u){
    double result = 0.0;
    result = u->x*u->x + u->y*u->y + u->z*u->z + u->w*u->w;
    result = sqrt(result+DBL_EPSILON);
    return result;
};

__device__
void la_d4_blas_dscal(double alpha, double4 *u){
    u->x = alpha*u->x;
    u->y = alpha*u->y;
    u->z = alpha*u->z;
    u->w = alpha*u->w;
};

__device__
void la_d4_blas_dscal(double alpha, double4 *uu, const size_t &col_size){
    for ( size_t row = 0; row < col_size; ++row ) {
        uu[row].x = alpha * uu[row].x;
        uu[row].y = alpha * uu[row].y;
        uu[row].z = alpha * uu[row].z;
        uu[row].w = alpha * uu[row].w;
    }
};

__device__
void la_d4_blas_dscal(const int i, double alpha, double4 *uu){
    uu[i].x = alpha*uu[i].x;
    uu[i].y = alpha*uu[i].y;
    uu[i].z = alpha*uu[i].z;
    uu[i].w = alpha*uu[i].w;
};

__device__
void la_d4_set(int i, double4 *uu,const double4 *v){
    uu[i].x = v->x;
    uu[i].y = v->y;
    uu[i].z = v->z;
    uu[i].w = v->w;
};
//#####################################################################################
__device__
void la_d4_add( double4 *u,const double4 *v){
    u->x = u->x + v->x;
    u->y = u->y + v->y;
    u->z = u->z + v->z;
    u->w = u->w + v->w;
};

__device__
void la_d4_add(int i, double4 *uu,const double4 *v){
    uu[i].x = uu[i].x + v->x;
    uu[i].y = uu[i].y + v->y;
    uu[i].z = uu[i].z + v->z;
    uu[i].w = uu[i].w + v->w;
};

__device__
void la_d4_add(double4 *uu, const double4 *vv, const size_t &col_size){
    for ( size_t row = 0; row < col_size; ++row ) {
        uu[row].x = uu[row].x + vv[row].x;
        uu[row].y = uu[row].y + vv[row].y;
        uu[row].z = uu[row].z + vv[row].z;
        uu[row].w = uu[row].w + vv[row].w;
    }
};
//#####################################################################################

/*
__device__
void la_d4_add(const int i, double4 *uu,const double4 *ww){
    uu[i].x = uu[i].x + ww[i].x;
    uu[i].y = uu[i].y + ww[i].y;
    uu[i].z = uu[i].z + ww[i].z;
    uu[i].w = uu[i].w + ww[i].w;
};
*/

__device__
void la_d4_sub( double4 *u,const double4 *v){
    u->x = u->x - v->x;
    u->y = u->y - v->y;
    u->z = u->z - v->z;
    u->w = u->w - v->w;
};

__device__
double la_d_pow_2(const double w){
    return w*w;
};

__device__
void la_d4_memcpy(double4 *u,const double4 *v){
    u->x = v->x;
    u->y = v->y;
    u->z = v->z;
    u->w = v->w;
};

__device__
const double4 la_d4_matrix_row(double4 *uu,const size_t &row_num) {
    const double4 v={uu[row_num].x, uu[row_num].y, uu[row_num].z, uu[row_num].w};
    return v;
};

__device__
double la_d4_vector_get(double4 *uu,const size_t &col_num){
    return uu[col_num].x;
};

__device__
double la_d_vector_get(const double *u,const size_t &col_num){
    return u[col_num];
};

__device__
void la_d4_sub(double4 *uu, const double4 *vv, const size_t &col_size){
    for ( size_t row = 0; row < col_size; ++row ) {
        uu[row].x = uu[row].x - vv[row].x;
        uu[row].y = uu[row].y - vv[row].y;
        uu[row].z = uu[row].z - vv[row].z;
        uu[row].w = uu[row].w - vv[row].w;
    }
};


__host__ __device__
void la_d4dealloc(double4 *uu){
    free((double*)uu);
};

__host__ __device__
void la_d4print(const double4 *uu,const size_t &col_size){
    for (size_t row = 0; row < col_size; ++row) {
        printf("\t%f\t%f\t%f\t%f", uu[row].x, uu[row].y, uu[row].z, uu[row].w);
        //printf("\t%f\t%f\t%f", uu[row].x, uu[row].y, uu[row].z);
        printf("\n");
    }
};

__host__ __device__
void la_d4print(const double4 *uu){
        printf("\t%f\t%f\t%f\t%f \n", uu->x, uu->y, uu->z, uu->w);
};
__host__ __device__
void la_d4print(const double4 u){
    printf("\t%f\t%f\t%f\t%f \n", u.x, u.y, u.z, u.w);
};


#endif
