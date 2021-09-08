//
// Created by ge56beh on 07.09.21.
//

#ifndef LA_UTILS_CU_H
#define LA_UTILS_CU_H

inline void la_arr_alloc( double **uu, const size_t &row_size, const size_t &col_size){
    uu = new double* [row_size];
    for ( size_t row = 0; row < row_size; ++row ) {
        uu[row] = new double[col_size];
    }
};

inline double** la_arr_calloc( const size_t &row_size, const size_t &col_size){
    double** uu = new double* [row_size];
    for ( size_t row = 0; row < row_size; ++row ) {
        uu[row] = new double[col_size];
        for (size_t col = 0; col < col_size; ++col) {
            uu[row][col] = 0.0;
        }
        return uu;
    }
};

inline void la_gslalloc( double **uu, const gsl_matrix *vv,  const size_t &row_size, const size_t &col_size){
    for ( size_t row = 0; row < row_size; ++row ) {
        for (size_t col = 0; col < col_size; ++col) {
            uu[row][col] = gsl_matrix_get(vv, row, col);
        }
    }
};

inline void la_dealloc( double **uu,  size_t &row_size){
    for (size_t row = 0; row < row_size; ++row) {
        delete[] uu[row];
    }
    delete[] uu;
};

inline void la_gsl_cross(const gsl_vector *a, const gsl_vector *b, gsl_vector *c){
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

inline void la_gsl_blas_dscal(double alpha, gsl_vector *x){

};

inline double* la_gsl_matrix_const_row(double **uu,const size_t &row_num,const size_t &row_size){
    double *v = new double [row_size];
    for(size_t col = 0; col<row_size; ++col){
        v[col] = uu[row_num][col];
    }
    return v;
};

inline double la_gsl_vector_get(double *u,const size_t &row_num){
    return u[row_num];
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

inline double la_gsl_pow_2(const double w){
    return w*w;
};

#endif
