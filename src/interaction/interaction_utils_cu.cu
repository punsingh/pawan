//
// Created by ge56beh on 13.09.21.
//

/*
#ifdef __CUDACC__
#define CUDA_HOSTDEV __host__ __device__
#endif
*/
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <cuda.h>
#include <cooperative_groups.h>
namespace cg = cooperative_groups;

#include "interaction_utils_cu.h"
#include "test_cuda.h"
#include "la_utils_cu.h"
#include "cuda_utils_cu.h"
#include "src/utils/gsl_utils.h"

#define CUDART_PI       3.14159265358979323846	/* pi */
#define CUDART_1_PI     0.31830988618379067154  /* 1/pi */
#define CUDART_SQRT1_2	0.70710678118654752440	/* 1/sqrt(2) */
#define BLOCK_SIZE      64                      /*max value shared memory-limited right now to 64*/

__device__
void cuKERNEL(	const double &rho,
                const double &sigma,
                double &q,
                double &F,
                double &Z){
    double rho_bar = rho/sigma;
    double sig3 = sigma*sigma*sigma;
    double phi = 0.25*CUDART_1_PI*erf(CUDART_SQRT1_2*rho_bar)/sig3;
    Z = 0.5*exp(-0.5*rho_bar*rho_bar)/sig3/pow(CUDART_PI,1.5);
    q = (phi/rho_bar - Z)/la_d_pow_2(rho_bar);
    F = (Z - 3*q)/la_d_pow_2(rho);
};

__device__
void cuVELOCITY(  const double &kernel,
                  const double4 *vorticity,
                  const double4 *displacement,
                  double4 *velocity){
    la_d4_cross(vorticity,displacement,velocity);
    la_d4_blas_dscal(kernel,velocity);
};

__device__
void cuVORSTRETCH( const double &q,
                   const double &F,
                   const double4 *source_vorticity,
                   const double4 *target_vorticity,
                   const double4 *displacement,
                   double4 *retvorcity){
    //trgXsrc = a_target x a_source
    double4 trgXsrc = {0.0, 0.0, 0.0, 0.0};
    la_d4_cross(target_vorticity,source_vorticity,&trgXsrc);

    //temp measure, needs fixing!!!
    if (target_vorticity->x == source_vorticity->x  &&  target_vorticity->y == source_vorticity->y  &&
        target_vorticity->z == source_vorticity->z  &&  target_vorticity->w == source_vorticity->w){
        trgXsrc = {0.0,0.0,0.0,0.0,};
    }

    // da/dt = q*(a_trg x a_src)
    double4 crossed = {q*trgXsrc.x, q*trgXsrc.y, q*trgXsrc.z, q*trgXsrc.w};

    // da/dt = F*[disp.(a_trg x a_src)]disp
    double roaxa = 0.0;
    la_d4_blas_ddot(displacement,&trgXsrc,&roaxa);
    double4 stretch = {displacement->x, displacement->y, displacement->z, displacement->w};
    la_d4_blas_dscal(F*roaxa,&stretch);

    la_d4_add(retvorcity,&crossed);
    la_d4_add(retvorcity,&stretch);
};

__device__
void cuDIFFUSION(	 const double nu,
                     const double &sigma,
                     const double &Z,
                     const double4 *source_vorticity,
                     const double4 *target_vorticity,
                     const double &source_volume,
                     const double &target_volume,
                     double4 *retvorcity){
    // va12 = volume_target*vorticity_source
    double4 va12 = {source_vorticity->x, source_vorticity->y, source_vorticity->z, source_vorticity->w};
    la_d4_blas_dscal(target_volume,&va12);

    // va21 = volume_source*vorticity_target
    double4 va21 = {target_vorticity->x, target_vorticity->y, target_vorticity->z, target_vorticity->w};
    la_d4_blas_dscal(source_volume,&va21);

    // dva = 2*nu*Z*(va12 - va21)/sigma^2
    double sig12 = 0.5*sigma*sigma;
    double4 dva = {va12.x, va12.y, va12.z, va12.w};
    la_d4_sub(&dva,&va21);
    la_d4_blas_dscal(Z*nu/sig12,&dva);

    // da = da + dva
    la_d4_add(retvorcity,&dva);
};

__device__
void cuINTERACT( const double nu,
                 const double &s_source,
                 const double &s_target,
                 const double4 *r_source,
                 const double4 *r_target,
                 const double4 *a_source,
                 const double4 *a_target,
                 const double &v_source,
                 const double &v_target,
                 double4 *dr_target,
                 double4 *da_target){
    // Kernel Computation
    double4 displacement = {r_target->x - r_source->x, r_target->y - r_source->y,
                            r_target->z - r_source->z,r_target->w - r_source->w};
    double rho = la_d4_blas_dnrm2_soft(&displacement);
    double q = 0.0, F = 0.0, Z = 0.0;
    double sigma = sqrt(s_source*s_source + s_target*s_target)/2.0;

    // Velocity computation
    double4 *dr = la_d4calloc(1);
    cuKERNEL(rho,sigma,q,F,Z);
    cuVELOCITY(-q,a_target,&displacement,dr);
    la_d4_add(dr_target,dr);

    // Rate of change of vorticity computation
    double4 *da =  la_d4calloc(1);
    cuVORSTRETCH(q,F,a_source,a_target,&displacement,da);
    cuDIFFUSION(nu,sigma,Z,a_source,a_target,v_source,v_target,da);
    la_d4_sub(da_target,da);

    la_d4dealloc(da);
    la_d4dealloc(dr);
}

__device__
void cuinteract(const int threadnum,
                const double nu,
                const double4 *Wpos_d4_arr,
                const double4 *Wvor_d4_arr,
                const double *Wrad_d_vec,
                const double *Wvol_d_vec,
                const int n,
                double4 *rates){

    const double4 r_src = Wpos_d4_arr[threadnum];
    const double4 a_src = Wvor_d4_arr[threadnum];
    const double s_src = Wrad_d_vec[threadnum];
    const double v_src = Wvol_d_vec[threadnum];

    for (int blocknum = 0; blocknum < gridDim.x; blocknum++) {

        __shared__ double4 wpos_d4_block[BLOCK_SIZE], wvor_d4_block[BLOCK_SIZE];
        __shared__ double wrad_d_block[BLOCK_SIZE], wvol_d_block[BLOCK_SIZE];

        wpos_d4_block[threadIdx.x] = Wpos_d4_arr[blocknum * blockDim.x + threadIdx.x];
        wvor_d4_block[threadIdx.x] = Wvor_d4_arr[blocknum * blockDim.x + threadIdx.x];
        wrad_d_block[threadIdx.x] = Wrad_d_vec[blocknum * blockDim.x + threadIdx.x];
        wvol_d_block[threadIdx.x] = Wvol_d_vec[blocknum * blockDim.x + threadIdx.x];
        __syncthreads();
        //#pragma unroll
        for (int j = 0; j < BLOCK_SIZE; j++) {
            idx = {threadnum, blocknum * blockDim.x + threadIdx.x};
            const double4 r_trg = wpos_d4_block[j];
            const double4 a_trg = wvor_d4_block[j];
            const double s_trg = wrad_d_block[j];
            const double v_trg = wvol_d_block[j];

            cuINTERACT(nu, s_src, s_trg, &r_src, &r_trg, &a_src, &a_trg, v_src, v_trg, rates, rates + 1);
        }
        __syncthreads();  //necessary before operations move onto next block
    }
}
/*!
 *  (Cuda kernel) Computes the time-integrated position and vorticity
 */
__global__
void cuINTERACT_rk4(double nu,
                    double4 *Wpos_d4_arr,
                    double4 *Wvor_d4_arr,
                    double4 *Wvel_d4_arr,
                    double4 *Wretvor_d4_arr,
                    double *Wrad_d_vec,
                    double *Wvol_d_vec,
                    double4 *Wpos_d4_arr_x1,
                    double4 *Wvor_d4_arr_x1,
                    double4 *Wpos_d4_arr_x2,
                    double4 *Wvor_d4_arr_x2,
                    double4 *Wpos_d4_arr_x3,
                    double4 *Wvor_d4_arr_x3,
                    const int n,
                    const double dt,
                    const size_t t_steps){
    cg::grid_group grid = cg::this_grid();
    //cg::thread_block block = cg::this_thread_block();

    int threadnum = blockDim.x * blockIdx.x + threadIdx.x;

    if (threadnum < n) {
        for (size_t t_n = 1; t_n<=t_steps; ++t_n) {

            double4 *k1 = la_d4calloc(2);
            double4 *k2 = la_d4calloc(2);
            double4 *k3 = la_d4calloc(2);
            double4 *k4 = la_d4calloc(2);

//######################################################################################################################
            cuinteract(threadnum, nu, Wpos_d4_arr, Wvor_d4_arr, Wrad_d_vec, Wvol_d_vec, n, k1);

            cg::sync(grid);
            la_d4_set(threadnum, Wpos_d4_arr_x1, k1);    la_d4_set(threadnum, Wvor_d4_arr_x1, k1+1);
            la_d4_blas_dscal(threadnum,0.5*dt, Wpos_d4_arr_x1);   la_d4_blas_dscal(threadnum,0.5*dt, Wvor_d4_arr_x1);
            la_d4_add(threadnum, Wpos_d4_arr_x1, Wpos_d4_arr+threadnum);    la_d4_add(threadnum, Wvor_d4_arr_x1, Wvor_d4_arr+threadnum);

            cg::sync(grid);
//######################################################################################################################
            cuinteract(threadnum, nu, Wpos_d4_arr_x1, Wvor_d4_arr_x1, Wrad_d_vec, Wvol_d_vec, n, k2);

            cg::sync(grid);
            la_d4_set(threadnum, Wpos_d4_arr_x2, k2);    la_d4_set(threadnum, Wvor_d4_arr_x2, k2+1);
            la_d4_blas_dscal(threadnum,0.5*dt, Wpos_d4_arr_x2);   la_d4_blas_dscal(threadnum,0.5*dt, Wvor_d4_arr_x2);
            la_d4_add(threadnum, Wpos_d4_arr_x2, Wpos_d4_arr+threadnum);    la_d4_add(threadnum, Wvor_d4_arr_x2, Wvor_d4_arr+threadnum);

            cg::sync(grid);
//######################################################################################################################
            cuinteract(threadnum, nu, Wpos_d4_arr_x2, Wvor_d4_arr_x2, Wrad_d_vec, Wvol_d_vec, n, k3);

            cg::sync(grid);
            la_d4_set(threadnum, Wpos_d4_arr_x3, k3);    la_d4_set(threadnum, Wvor_d4_arr_x3, k3+1);
            la_d4_blas_dscal(threadnum,dt, Wpos_d4_arr_x3);   la_d4_blas_dscal(threadnum,dt, Wvor_d4_arr_x3);
            la_d4_add(threadnum, Wpos_d4_arr_x3, Wpos_d4_arr+threadnum);    la_d4_add(threadnum, Wvor_d4_arr_x3, Wvor_d4_arr+threadnum);

            cg::sync(grid);
//######################################################################################################################
            cuinteract(threadnum, nu, Wpos_d4_arr_x3, Wvor_d4_arr_x3, Wrad_d_vec, Wvol_d_vec, n, k4);

            cg::sync(grid);
            la_d4_add(k1,k4,2);   la_d4_blas_dscal(dt/6.0,k1,2);
            la_d4_add(k2,k3,2);   la_d4_blas_dscal(dt/3.0,k2,2);
            la_d4_add(k1,k2,2);
            la_d4_add(threadnum, Wpos_d4_arr, k1);    la_d4_add(threadnum, Wvor_d4_arr, k1+1);
            la_d4_set(threadnum, Wvel_d4_arr, k4);    la_d4_set(threadnum, Wretvor_d4_arr, k4 + 1);

            cg::sync(grid);
//######################################################################################################################

            la_d4dealloc(k1); la_d4dealloc(k2); la_d4dealloc(k3); la_d4dealloc(k4);

            if (threadnum == 1){ printf("Step %d \n", t_n);}
        }
   }
};

void cuda_main(double &nu,
               double **Wpos_arr,
               double **Wvor_arr,
               double **Wvel_arr,
               double **Wretvor_arr,
               double *Wrad_vec,
               double *Wvol_vec,
               const size_t numParticles,
               const size_t numDimensions,
               const double dt,
               const size_t t_steps){

    int num_gpu = 0;  // number of CUDA GPUs
    cudaGetDeviceCount(&num_gpu); //get number of gpus available
    if (num_gpu < 1) {
        printf("no CUDA capable GPUs were detected \n");
        return;
    } else {
        printf("%d CUDA capable GPUs were detected \n", num_gpu);
    }
    int gpuID = num_gpu - 1; //the last (non-default) gpu is 'usually' free
    if (cudaSetDevice(gpuID) != cudaSuccess)
        printf("something went wrong setting gpu num %d \n", gpuID);
    cudaDeviceProp cudprop;
    cudaDeviceProperties_print(&cudprop, gpuID);
    printf(" \n \n ");

    int d4_bytes = numParticles * sizeof(double4);    //double4 type heap allocation on device
    int d_bytes = numParticles * sizeof(double);     //double type heap allocation on device

    double4 *Wpos_d4_arr = la_to_d4alloc(Wpos_arr, numParticles, numDimensions);
    double4 *dev_wpos;
    cudaMalloc(&dev_wpos, d4_bytes);
    cudaMemcpy(dev_wpos, Wpos_d4_arr, d4_bytes, cudaMemcpyHostToDevice);

    double4 *Wvor_d4_arr = la_to_d4alloc(Wvor_arr, numParticles, numDimensions);
    double4 *dev_wvor;
    cudaMalloc(&dev_wvor, d4_bytes);
    cudaMemcpy(dev_wvor, Wvor_d4_arr, d4_bytes, cudaMemcpyHostToDevice);

    double4 *Wvel_d4_arr = la_to_d4alloc(Wvel_arr, numParticles, numDimensions);
    double4 *dev_wvel;
    cudaMalloc(&dev_wvel, d4_bytes);
    cudaMemcpy(dev_wvel, Wvel_d4_arr, d4_bytes, cudaMemcpyHostToDevice);

    double4 *Wretvor_d4_arr = la_to_d4alloc(Wretvor_arr, numParticles, numDimensions);
    double4 *dev_wretvor;
    cudaMalloc(&dev_wretvor, d4_bytes);
    cudaMemcpy(dev_wretvor, Wretvor_d4_arr, d4_bytes, cudaMemcpyHostToDevice);

    double4 *dev_wpos_x1, *dev_wpos_x2, *dev_wpos_x3, *dev_wpos_x4;
    cudaMalloc(&dev_wpos_x1, d4_bytes); cudaMalloc(&dev_wpos_x2, d4_bytes); cudaMalloc(&dev_wpos_x3, d4_bytes); cudaMalloc(&dev_wpos_x4, d4_bytes);
    cudaMemset(dev_wpos_x1,0.0,d4_bytes); cudaMemset(dev_wpos_x2,0.0,d4_bytes); cudaMemset(dev_wpos_x3,0.0,d4_bytes); cudaMemset(dev_wpos_x4,0.0,d4_bytes);

    double4 *dev_wvor_x1, *dev_wvor_x2, *dev_wvor_x3, *dev_wvor_x4;
    cudaMalloc(&dev_wvor_x1, d4_bytes); cudaMalloc(&dev_wvor_x2, d4_bytes); cudaMalloc(&dev_wvor_x3, d4_bytes); cudaMalloc(&dev_wvor_x4, d4_bytes);
    cudaMemset(dev_wvor_x1,0.0,d4_bytes); cudaMemset(dev_wvor_x2,0.0,d4_bytes); cudaMemset(dev_wvor_x3,0.0,d4_bytes); cudaMemset(dev_wvor_x4,0.0,d4_bytes);

    double *dev_wrad;
    cudaMalloc(&dev_wrad, d_bytes);
    cudaMemcpy(dev_wrad, Wrad_vec, d_bytes, cudaMemcpyHostToDevice);

    double *dev_wvol;
    cudaMalloc(&dev_wvol, d_bytes);
    cudaMemcpy(dev_wvol, Wvol_vec, d_bytes, cudaMemcpyHostToDevice);

    int nBlocks = (numParticles + BLOCK_SIZE - 1) / BLOCK_SIZE;

    //cuda kernel arguments
    void *kernelArgs[] = {
            (void *)&nu,  (void *)&dev_wpos, (void *)&dev_wvor, (void *)&dev_wvel,
            (void *)&dev_wretvor, (void *)&dev_wrad, (void *)&dev_wvol,
            (void *)&dev_wpos_x1, (void *)&dev_wvor_x1,
            (void *)&dev_wpos_x2, (void *)&dev_wvor_x2,  (void *)&dev_wpos_x3, (void *)&dev_wvor_x3,
            (void *)&numParticles, (void *)&dt,  (void *)&t_steps
    };
    dim3 dimGrid(nBlocks,1,1);  //check cuda-samples (multiple blocks possible per SM)
    dim3 dimBlock(BLOCK_SIZE,1,1);

    //lauching cuda kernels with 'cooperative groups'
    cudaLaunchCooperativeKernel((void *) cuINTERACT_rk4, dimGrid, dimBlock, kernelArgs);

    //cudaDeviceSynchronize();
    cudaMemcpy(Wpos_d4_arr, dev_wpos, d4_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(Wvor_d4_arr, dev_wvor, d4_bytes, cudaMemcpyDeviceToHost);
    //cudaMemcpy(Wrad_vec, dev_wrad, d4_bytes, cudaMemcpyDeviceToHost);
    //cudaMemcpy(Wvol_vec, dev_wvol, d4_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(Wvel_d4_arr, dev_wvel, d4_bytes, cudaMemcpyDeviceToHost);
    cudaMemcpy(Wretvor_d4_arr, dev_wretvor, d4_bytes, cudaMemcpyDeviceToHost);

    printf("After CUDA---------------------- \n");
    printf("Wpos_arr-----------------\n");
    la_d4print(Wpos_d4_arr, numParticles);
    printf("Wvor_arr-----------------\n");
    la_d4print(Wvor_d4_arr, numParticles);
    printf("Wvel_arr-----------------\n");
    la_d4print(Wvel_d4_arr, numParticles);
    printf("Wretvor_arr-----------------\n");
    la_d4print(Wretvor_d4_arr, numParticles);

    la_to_d4alloc_reverse(Wpos_arr,Wpos_d4_arr ,numParticles,numDimensions);
    la_to_d4alloc_reverse(Wvor_arr,Wvor_d4_arr ,numParticles,numDimensions);
    la_to_d4alloc_reverse(Wvel_arr,Wvel_d4_arr ,numParticles,numDimensions);
    la_to_d4alloc_reverse(Wretvor_arr,Wretvor_d4_arr ,numParticles,numDimensions);

    free(Wpos_d4_arr);
    free(Wvor_d4_arr);
    free(Wvel_d4_arr);
    free(Wretvor_d4_arr);

    cudaFree(dev_wpos);
    cudaFree(dev_wvor);
    cudaFree(dev_wvel);
    cudaFree(dev_wretvor);
    cudaFree(dev_wrad);
    cudaFree(dev_wvol);
}

