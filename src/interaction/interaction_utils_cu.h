//
// Created by Sumeet Kumar on 26/08/21.
//

#ifndef PAWAN_INTERACTION_UTILS_CU_H
#define PAWAN_INTERACTION_UTILS_CU_H



__device__ void KERNEL(	const double &rho,
                           const double &sigma,
                           double &q,
                           double &F,
                           double &Z){
    double rho_bar = rho/sigma;
    double sig3 = sigma*sigma*sigma;
    double phi = 0.25*M_1_PI*erf(M_SQRT1_2*rho_bar)/sig3;
    Z = 0.5*exp(-0.5*rho_bar*rho_bar)/sig3/pow(M_PI,1.5);
    q = (phi/rho_bar - Z)/gsl_pow_2(rho_bar);
    F = (Z - 3*q)/gsl_pow_2(rho);
};

__device__ void VELOCITY(	const double &kernel,
                             const gsl_vector *vorticity,
                             const gsl_vector *displacement,
                             gsl_vector *velocity ){
    cugsl_cross(vorticity,displacement,velocity);
    cugsl_blas_dscal(kernel,velocity);
};

__device__ void INTERACT(	const double &nu,
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













#endif //PAWAN_INTERACTION_UTILS_CU_H
