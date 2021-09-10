//
// Created by Sumeet Kumar on 26/08/21.
//

#ifndef INTERACTION_UTILS_CU_H
#define INTERACTION_UTILS_CU_H

#define CUDART_PI       3.14159265358979323846	/* pi */
#define CUDART_1_PI     0.31830988618379067154  /* 1/pi */
#define CUDART_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */

#include <math.h>    //comment this when running CUDA
#include "la_utils_cu.h"
#include "math_utils_cu.h"
#include "src/utils/gsl_utils.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

/*
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
*/
//__host__ __device__
void KERNEL(	const double &rho,
                                    const double &sigma,
                                    double &q,
                                    double &F,
                                    double &Z){
    double rho_bar = rho/sigma;
    double sig3 = sigma*sigma*sigma;
    double phi = 0.25*CUDART_1_PI*erf(CUDART_SQRT1_2*rho_bar)/sig3;
    Z = 0.5*exp(-0.5*rho_bar*rho_bar)/sig3/pow(CUDART_PI,1.5);
    q = (phi/rho_bar - Z)/la_gsl_pow_2(rho_bar);
    F = (Z - 3*q)/la_gsl_pow_2(rho);
};

//__host__ __device__
void VELOCITY(	const double &kernel,

                                      const double *vorticity,
                                      const double *displacement,
                                     double *velocity
/*
                                      const gsl_vector *vorticity_gsl,
                                      const gsl_vector *displacement_gsl,
                                      gsl_vector *velocity_gsl
*/){
/*
    double* vorticity = la_gslalloc(vorticity_gsl,3);//----------------------
    double* displacement = la_gslalloc(displacement_gsl,3);//----------------------
    double* velocity = la_gslalloc(velocity_gsl,3);//----------------------
*/
    int ndim = 3;

    la_gsl_cross(vorticity,displacement,velocity,ndim);
    la_gsl_blas_dscal(kernel,velocity,ndim);
//    la_gslalloc_reverse(velocity,velocity_gsl,ndim);//-----------------
};

/*!
 * Compute rate of change of vorticity due to vorterx stretching
 * \param	q			q/rho kernel
 * \param	F			F kernel
 * \param	source_vorticity	gsl vector source vorticity
 * \param	target_vorticity	gsl vector target vorticity
 * \param	displacement		gsl vector displacement between source and target
 * \param	retvorcity		gsl vector output rate of change of vorticity
 */

//__host__ __device__
inline void VORSTRETCH(	const double &q,
                                               const double &F,

                                               const double *source_vorticity,
                                               const double *target_vorticity,
                                               const double *displacement,
                                               double *retvorcity
/*
                                               const gsl_vector *source_vorticity_gsl,
                                               const gsl_vector *target_vorticity_gsl,
                                               const gsl_vector *displacement_gsl,
                                               gsl_vector *retvorcity_gsl
*/                           ){

/*
    double* source_vorticity = la_gslalloc(source_vorticity_gsl,3);//----------------------
    double* target_vorticity = la_gslalloc(target_vorticity_gsl,3);//----------------------
    double* displacement = la_gslalloc(displacement_gsl,3);//----------------------
*/
    int ndim = 3;

    //trgXsrc = a_target x a_source
    double *trgXsrc = new double[ndim];
    la_gsl_cross(target_vorticity,source_vorticity,trgXsrc,ndim);

    // da/dt = q*(a_trg x a_src)
    double *crossed = new double[ndim];
    la_gsl_vector_memcpy(crossed, trgXsrc, ndim);
    la_gsl_blas_dscal(q,crossed,ndim);

    // da/dt = F*[disp.(a_trg x a_src)]disp
    double roaxa = 0.0;
    la_gsl_blas_ddot(displacement,trgXsrc,&roaxa, ndim);
    double *stretch = new double[ndim];
    la_gsl_vector_memcpy(stretch,displacement, ndim);
    la_gsl_blas_dscal(F*roaxa,stretch,ndim);
/*
    double *retvorcity = new double[ndim];//-----------------------------------
    for(size_t col = 0; col<ndim; ++col){//-----------------------------------
        retvorcity[col] = 0;//-----------------------------------
    }//-----------------------------------
*/
    la_gsl_vector_add(retvorcity,crossed,ndim);
    la_gsl_vector_add(retvorcity,stretch,ndim);
/*
    gsl_vector *retvorcity_gsl_tmp = gsl_vector_calloc(ndim);//-----------------------------------
    la_gslalloc_reverse(retvorcity,retvorcity_gsl_tmp,ndim);//-----------------
    gsl_vector_add(retvorcity_gsl,retvorcity_gsl_tmp);//-----------------
    delete[] retvorcity;//-----------------------------------------------
    gsl_vector_free(retvorcity_gsl_tmp);//-----------------------
*/
    delete[] trgXsrc;
    delete[] crossed;
    delete[] stretch;
};


/*! \fn inline void DIFFUSION(	const double &nu, const double &sigma, const double &Z, const gsl_vector *source_vorticity, const gsl_vector *target_vorticity, const double &source_volume, const double &target_volume, gsl_vector *retvorcity )
 * \brief Compute rate of change of vorticity due to viscous diffusion
 * \param	nu			double viscosity
 * \param	sigma			double smoothing radius
 * \param	Z			double Z kernel
 * \param	source_vorticity	gsl vector source vorticity
 * \param	target_vorticity	gsl vector target vorticity
 * \param	source_volume		double source volume
 * \param	target_volume		double target volume
 * \param	retvorcity		gsl vector output rate of change of vorticity
 */
//__host__ __device__
inline void DIFFUSION(	const double &nu,
                                              const double &sigma,
                                              const double &Z,
                                              const double *source_vorticity,
                                              const double *target_vorticity,
                                              //const gsl_vector *source_vorticity_gsl,
                                              //const gsl_vector *target_vorticity_gsl,
                                              const double &source_volume,
                                              const double &target_volume,
                                             double *retvorcity
                                             // gsl_vector *retvorcity
                                              ){
//    double* source_vorticity = la_gslalloc(source_vorticity_gsl,3);//----------------------
//    double* target_vorticity = la_gslalloc(target_vorticity_gsl,3);//----------------------

    int ndim = 3;

    // va12 = volume_target*vorticity_source
    double *va12 = new double[ndim];
    la_gsl_vector_memcpy(va12,source_vorticity,ndim);
    la_gsl_blas_dscal(target_volume,va12,ndim);

    // va21 = volume_source*vorticity_target
    double *va21 = new double[ndim];
    la_gsl_vector_memcpy(va21,target_vorticity,ndim);
    la_gsl_blas_dscal(source_volume,va21,ndim);

    // dva = 2*nu*Z*(va12 - va21)/sigma^2
    double sig12 = 0.5*sigma*sigma;
    double *dva = new double[ndim];
    la_gsl_vector_memcpy(dva, va12, ndim);
    la_gsl_vector_sub(dva,va21,ndim);
    la_gsl_blas_dscal(Z*nu/sig12,dva,ndim);

    // da = da + dva
    la_gsl_vector_add(retvorcity,dva,3);
 /*
    gsl_vector *dva_gsl = gsl_vector_calloc(ndim);//---------------------------
    la_gslalloc_reverse(dva, dva_gsl, ndim);//---------------------------
    gsl_vector_add(retvorcity,dva_gsl);//---------------------------
    gsl_vector_free(dva_gsl);//---------------------------
 */
    delete[] va21;
    delete[] va12;
    delete[] dva;

};

//__host__


inline void INTERACT_(	const double &nu,
                         const double &s_source,
                         const double &s_target,
                         const double *r_source,
                         const double *r_target,
                         const double *a_source,
                         const double *a_target,
                         const double &v_source,
                         const double &v_target,
                         double *dr_target,
                         double *da_target,
                         double &vx_source,
                         double &vy_source,
                         double &vz_source,
                         double &qx_source,
                         double &qy_source,
                         double &qz_source){
/*
    double* r_source = la_gslalloc(r_source_gsl,3);//----------------------
    double* r_target = la_gslalloc(r_target_gsl,3);//----------------------
    double* a_source = la_gslalloc(a_source_gsl,3);//----------------------
    double* a_target = la_gslalloc(a_target_gsl,3);//----------------------
*/
    // Kernel Computation
int ndim = 3;

double *displacement = new double[ndim];
la_gsl_vector_memcpy(displacement,r_target,ndim);
    la_gsl_vector_sub(displacement,r_source,ndim);
    double rho = la_gsl_blas_dnrm2(displacement,ndim);
    double q = 0.0, F = 0.0, Z = 0.0;
    double sigma = sqrt(la_gsl_pow_2(s_source) + la_gsl_pow_2(s_target))/2.0;
// Velocity computation
    double *dr = new double[ndim];
    for(size_t col = 0; col<ndim; ++col){
        dr[col] = 0;
    }
// Target
    KERNEL(rho,sigma,q,F,Z);
    VELOCITY(q,a_source,displacement,dr);
    la_gsl_vector_add(dr_target,dr, ndim);
/*
    gsl_vector *dr_gsl_tmp = gsl_vector_calloc(ndim);//-----------------------------------
    la_gslalloc_reverse(dr,dr_gsl_tmp,ndim);//-----------------
    gsl_vector_add(dr_target_gsl,dr_gsl_tmp);//-----------------
    gsl_vector_free(dr_gsl_tmp);//-----------------------
*/
// Source
    VELOCITY(-q,a_target,displacement,dr);
    vx_source = dr[0];
    vy_source = dr[1];
    vz_source = dr[2];

    double *da = new double [ndim];
    for(size_t col = 0; col<ndim; ++col){
        da[col] = 0;
    }
    VORSTRETCH(q,F,a_source,a_target,displacement,da);
    DIFFUSION(nu,sigma,Z,a_source,a_target,v_source,v_target,da);
    la_gsl_vector_add(da_target,da,ndim);
// Target
// Source
    qx_source = -da[0];
    qy_source = -da[1];
    qz_source = -da[2];
/*
    gsl_vector *da_gsl = gsl_vector_calloc(ndim);//-----------------------------------
    la_gslalloc_reverse(da,da_gsl,ndim);//-----------------
    gsl_vector_add(da_target_gsl,da_gsl);//-----------------
    gsl_vector_free(da_gsl);//-----------------------
*/
    delete[] dr;
    delete[] da;
    delete[] displacement;

/*
    delete[] r_source;//-----------------------
    delete[] r_target;//-----------------------
    delete[] a_source;//-----------------------
    delete[] a_target;//-----------------------
*/


};


inline void INTERACT(	const double &nu,
                  const double &s_source,
                  const double &s_target,
                  const gsl_vector *r_source_gsl,
                  const gsl_vector *r_target_gsl,
                  const gsl_vector *a_source_gsl,
                  const gsl_vector *a_target_gsl,
                  const double &v_source,
                  const double &v_target,
                  gsl_vector *dr_target_gsl,
                  gsl_vector *da_target_gsl,
                  double &vx_source,
                  double &vy_source,
                  double &vz_source,
                  double &qx_source,
                  double &qy_source,
                  double &qz_source){

    double* r_source = la_gslalloc(r_source_gsl,3);//----------------------
    double* r_target = la_gslalloc(r_target_gsl,3);//----------------------
    double* a_source = la_gslalloc(a_source_gsl,3);//----------------------
    double* a_target = la_gslalloc(a_target_gsl,3);//----------------------

    // Kernel Computation
    int ndim = 3;

    double *displacement = new double[ndim];
    la_gsl_vector_memcpy(displacement,r_target,ndim);
    la_gsl_vector_sub(displacement,r_source,ndim);
    double rho = la_gsl_blas_dnrm2(displacement,ndim);
    double q = 0.0, F = 0.0, Z = 0.0;
    double sigma = sqrt(la_gsl_pow_2(s_source) + la_gsl_pow_2(s_target))/2.0;
// Velocity computation
    double *dr = new double[ndim];
    for(size_t col = 0; col<ndim; ++col){
        dr[col] = 0;
    }
// Target
    KERNEL(rho,sigma,q,F,Z);
    VELOCITY(q,a_source,displacement,dr);

    gsl_vector *dr_gsl_tmp = gsl_vector_calloc(ndim);//-----------------------------------
    la_gslalloc_reverse(dr,dr_gsl_tmp,ndim);//-----------------
    gsl_vector_add(dr_target_gsl,dr_gsl_tmp);//-----------------
    gsl_vector_free(dr_gsl_tmp);//-----------------------

// Source
    VELOCITY(-q,a_target,displacement,dr);
    vx_source = dr[0];
    vy_source = dr[1];
    vz_source = dr[2];

    double *da = new double [ndim];
    for(size_t col = 0; col<ndim; ++col){
        da[col] = 0;
    }
    VORSTRETCH(q,F,a_source,a_target,displacement,da);
    DIFFUSION(nu,sigma,Z,a_source,a_target,v_source,v_target,da);
// Target
// Source
    qx_source = -da[0];
    qy_source = -da[1];
    qz_source = -da[2];

    gsl_vector *da_gsl_tmp = gsl_vector_calloc(ndim);//-----------------------------------
    la_gslalloc_reverse(da,da_gsl_tmp,ndim);//-----------------
    gsl_vector_add(da_target_gsl,da_gsl_tmp);//-----------------
    gsl_vector_free(da_gsl_tmp);//-----------------------

    delete[] dr;
    delete[] da;
    delete[] displacement;


    delete[] r_source;//-----------------------
    delete[] r_target;//-----------------------
    delete[] a_source;//-----------------------
    delete[] a_target;//-----------------------



};



/*
__global__ void INTERACT(	const double &nu,
                             const double &s_source,
                             const double &s_target,
                             const gsl_vector *r_source,
                             const gsl_vector *r_target,
                             const gsl_vector *a_source,
                             const gsl_vector *a_target,
                             const double &v_source,
                             const double &v_target,
                             gsl_vector *dr_target,
                             gsl_vector *da_target,
                             double &vx_source,
                             double &vy_source,
                             double &vz_source,
                             double &qx_source,
                             double &qy_source,
                             double &qz_source){




    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) {
        float i_vx = 0.0f; float i_vy = 0.0f; float i_vz = 0.0f;
        float i_qx = 0.0f; float i_qy = 0.0f; float i_qz = 0.0f;

        for (int block_num = 0; block_num < gridDim.x; block_num++) {
            __shared__ float3 j_pos[BLOCK_SIZE];
            __shared__ float3 j_rad[BLOCK_SIZE];
            __shared__ float3 j_vol[BLOCK_SIZE];
            __shared__ float3 j_vor[BLOCK_SIZE];

            float4 tpos = p[block_num * blockDim.x + threadIdx.x];
            j_pos[threadIdx.x] = make_float3(tpos.x, tpos.y, tpos.z);
            __syncthreads();

#pragma unroll
            for (int j = 0; j < BLOCK_SIZE; j++) {
                float dx = spos[j].x - p[i].x;
                float dy = spos[j].y - p[i].y;
                float dz = spos[j].z - p[i].z;
                float distSqr = dx*dx + dy*dy + dz*dz + SOFTENING;
                float invDist = rsqrtf(distSqr);
                float invDist3 = invDist * invDist * invDist;

                Fx += dx * invDist3; Fy += dy * invDist3; Fz += dz * invDist3;
            }
            __syncthreads();
        }
        Wvel[i] += i_vx; Wvel[i] += i_vy; Wvel[i] += i_vz;
        Wretvor[i] += q_vx; Wretvor[i] += q_vy; Wretvor[i] += q_vz;
    }
}


// Kernel Computation
gsl_vector *displacement = gsl_vector_calloc(3);
gsl_vector_memcpy(displacement,r_target);
gsl_vector_sub(displacement,r_source);
double rho = gsl_blas_dnrm2(displacement);
double q = 0.0, F = 0.0, Z = 0.0;
double sigma = sqrt(gsl_pow_2(s_source) + gsl_pow_2(s_target))/2.0;

// Velocity computation
gsl_vector *dr = gsl_vector_calloc(3);
// Target
KERNEL(rho,sigma,q,F,Z);
VELOCITY(q,a_source,displacement,dr);
gsl_vector_add(dr_target,dr);
// Source
VELOCITY(-q,a_target,displacement,dr);
vx_source = gsl_vector_get(dr,0);
vy_source = gsl_vector_get(dr,1);
vz_source = gsl_vector_get(dr,2);

// Rate of change of vorticity computation
gsl_vector *da = gsl_vector_calloc(3);
VORSTRETCH(q,F,a_source,a_target,displacement,da);
DIFFUSION(nu,sigma,Z,a_source,a_target,v_source,v_target,da);
// Target
gsl_vector_add(da_target,da);
// Source
qx_source = -gsl_vector_get(da,0);
qy_source = -gsl_vector_get(da,1);
qz_source = -gsl_vector_get(da,2);

// Clean up
gsl_vector_free(dr);
gsl_vector_free(da);
gsl_vector_free(displacement);
};


*/










#endif //INTERACTION_UTILS_CU_H
