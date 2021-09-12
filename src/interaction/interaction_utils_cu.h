//
// Created by Sumeet Kumar on 26/08/21.
//

#ifndef INTERACTION_UTILS_CU_H
#define INTERACTION_UTILS_CU_H

#include <math.h>    //comment this when running CUDA
#include "la_utils_cu.h"
#include "math_utils_cu.h"
//#include "cuda_utils_cu.h"
#include "src/utils/gsl_utils.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

#define CUDART_PI       3.14159265358979323846	/* pi */
#define CUDART_1_PI     0.31830988618379067154  /* 1/pi */
#define CUDART_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */


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
){

    int ndim = 3;

    la_gsl_cross(vorticity,displacement,velocity,ndim);
    la_gsl_blas_dscal(kernel,velocity,ndim);
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
                           ){

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

    la_gsl_vector_add(retvorcity,crossed,ndim);
    la_gsl_vector_add(retvorcity,stretch,ndim);

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
                                              const double &source_volume,
                                              const double &target_volume,
                                             double *retvorcity
                                              ){

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

    delete[] va21;
    delete[] va12;
    delete[] dva;
};



//__global__
inline void INTERACT(	const double &nu,
                         const double &s_source,
                         const double &s_target,
                         const double *r_source,
                         const double *r_target,
                         const double *a_source,
                         const double *a_target,
                         const double &v_source,
                         const double &v_target,
                         double *dr_target,
                         double *da_target){
    size_t ndim = 3;

    // Kernel Computation
    double *displacement = new double[ndim];
    la_gsl_vector_memcpy(displacement,r_target,ndim);
    la_gsl_vector_sub(displacement,r_source,ndim);
    double rho = la_gsl_blas_dnrm2_soft(displacement,ndim);
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
    // Source
//    VELOCITY(-q,a_target,displacement,dr);

    // Rate of change of vorticity computation
    double *da = new double [ndim];
    for(size_t col = 0; col<ndim; ++col){
        da[col] = 0;
    }
    VORSTRETCH(q,F,a_source,a_target,displacement,da);
    DIFFUSION(nu,sigma,Z,a_source,a_target,v_source,v_target,da);
    // Target
    la_gsl_vector_add(da_target,da,ndim);

    delete[] dr;
    delete[] da;
    delete[] displacement;
};

/*
//__host__
inline void INTERACT(	const double &nu,
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
    int ndim = 3;

    // Kernel Computation
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
    // Source
    VELOCITY(-q,a_target,displacement,dr);
    vx_source = dr[0];
    vy_source = dr[1];
    vz_source = dr[2];

    // Rate of change of vorticity computation
    double *da = new double [ndim];
    for(size_t col = 0; col<ndim; ++col){
        da[col] = 0;
    }
    VORSTRETCH(q,F,a_source,a_target,displacement,da);
    DIFFUSION(nu,sigma,Z,a_source,a_target,v_source,v_target,da);
    // Target
    la_gsl_vector_add(da_target,da,ndim);
    // Source
    qx_source = -da[0];
    qy_source = -da[1];
    qz_source = -da[2];

    delete[] dr;
    delete[] da;
    delete[] displacement;
};
*/

#endif //INTERACTION_UTILS_CU_H
