//
// Created by Sumeet Kumar on 26/08/21.
//

#ifndef INTERACTION_UTILS_CU_H
#define INTERACTION_UTILS_CU_H

#define CUDART_PI_F 3.141592654f

#include <math.h>    //comment this when running CUDA
#include "la_utils_cu.h"
#include "math_utils_cu.h"
#include "src/utils/gsl_utils.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>

//__host__ __device__
void KERNEL(	const double &rho,
                                    const double &sigma,
                                    double &q,
                                    double &F,
                                    double &Z){
    double rho_bar = rho/sigma;
    double sig3 = sigma*sigma*sigma;
    double phi = 0.25*(1/CUDART_PI_F)*erf(sqrt(0.5)*rho_bar)/sig3;
    Z = 0.5*exp(-0.5*rho_bar*rho_bar)/sig3/pow(CUDART_PI_F,1.5);
    q = (phi/rho_bar - Z)/la_gsl_pow_2(rho_bar);
    F = (Z - 3*q)/la_gsl_pow_2(rho);
};

//__host__ __device__
void VELOCITY(	const double &kernel,
                                      const gsl_vector *vorticity,
                                      const gsl_vector *displacement,
                                      gsl_vector *velocity ){
    gsl_cross(vorticity,displacement,velocity);
    gsl_blas_dscal(kernel,velocity);
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
                                               const gsl_vector *source_vorticity,
                                               const gsl_vector *target_vorticity,
                                               const gsl_vector *displacement,
                                               gsl_vector *retvorcity ){
    // a_target x a_source
    gsl_vector *trgXsrc = gsl_vector_calloc(3);
    gsl_cross(target_vorticity,source_vorticity,trgXsrc);

    gsl_vector *crossed = gsl_vector_calloc(3);
    gsl_vector *stretch = gsl_vector_calloc(3);

    // da/dt = q*(a_trg x a_src)
    gsl_vector_memcpy(crossed,trgXsrc);
    gsl_blas_dscal(q,crossed);

    // da/dt = F*[disp.(a_trg x a_src)]disp
    double roaxa = 0.0;;
    gsl_blas_ddot(displacement,trgXsrc,&roaxa);
    gsl_vector_memcpy(stretch,displacement);
    gsl_blas_dscal(F*roaxa,stretch);

    //gsl_vector_set_zero(retvorcity);
    gsl_vector_add(retvorcity,crossed);
    gsl_vector_add(retvorcity,stretch);

    gsl_vector_free(trgXsrc);
    gsl_vector_free(crossed);
    gsl_vector_free(stretch);

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
                                              const gsl_vector *source_vorticity,
                                              const gsl_vector *target_vorticity,
                                              const double &source_volume,
                                              const double &target_volume,
                                              gsl_vector *retvorcity ){

    gsl_vector *va12 = gsl_vector_calloc(3);
    gsl_vector *va21 = gsl_vector_calloc(3);
    gsl_vector *dva = gsl_vector_calloc(3);

    // va12 = volume_target*vorticity_source
    gsl_vector_memcpy(va12,source_vorticity);
    gsl_blas_dscal(target_volume,va12);

    // va21 = volume_source*vorticity_target
    gsl_vector_memcpy(va21,target_vorticity);
    gsl_blas_dscal(source_volume,va21);

    // dva = 2*nu*Z*(va12 - va21)/sigma^2
    double sig12 = 0.5*sigma*sigma;
    gsl_vector_memcpy(dva,va12);
    gsl_vector_sub(dva,va21);
    gsl_blas_dscal(Z*nu/sig12,dva);

    // da = da + dva
    gsl_vector_add(retvorcity,dva);

    gsl_vector_free(va12);
    gsl_vector_free(va21);
    gsl_vector_free(dva);
};

//__host__
inline void INTERACT(	const double &nu,
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
//double *displacement = new double[3];
//la_gsl_vector_memcpy(displacement,r_target,W->_numDimensions);
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
