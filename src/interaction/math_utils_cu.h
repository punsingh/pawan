//
// Created by Sumeet Kumar on 26/08/21.
//

#ifndef PAWAN_MATH_UTILS_CU_H
#define PAWAN_MATH_UTILS_CU_H



__device__ double erf ( double  x ){}

__device__ double pow ( double  x, double  y ){}

__device__ double sqrt ( double  x ){}

__device__ void cugsl_cross(const gsl_vector *a, const gsl_vector *b, gsl_vector *c){
    if(a->size==3 && b->size==3 && c->size==3){
        double ax = gsl_vector_get(a,0);
        double ay = gsl_vector_get(a,1);
        double az = gsl_vector_get(a,2);
        double bx = gsl_vector_get(b,0);
        double by = gsl_vector_get(b,1);
        double bz = gsl_vector_get(b,2);
        gsl_vector_set(c,0,ay*bz-az*by);
        gsl_vector_set(c,1,az*bx-ax*bz);
        gsl_vector_set(c,2,ax*by-ay*bx);
    }
};

__device__ void cugsl_blas_dscal(double *scale,velocity){}


#endif //PAWAN_MATH_UTILS_CU_H
