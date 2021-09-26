//
// Created by Sumeet Kumar on 26/08/21.
//

#ifndef INTERACTION_UTILS_CU_H
#define INTERACTION_UTILS_CU_H

#include "interaction_utils_cu.h"
#include <math.h>    //comment this when running CUDA
#include "la_utils_cu.h"
#include "cuda_utils_cu.h"
#include "src/utils/gsl_utils.h"
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <stdio.h>
#include <cuda.h>

void cuda_main(   double &nu,
               double **Wpos_arr,
               double **Wvor_arr,
               double** Wvel_arr,
               double** Wretvor_arr,
               double* Wrad_vec,
               double* Wvol_vec,
               const size_t numParticles,
               const size_t numDimensions,
               const double dt,
               const size_t t_steps);
/*
void cuda_main(double &nu,
                   double **Wpos_arr,
                   double **Wvor_arr,
                   double** Wvel_arr,
                   double** Wretvor_arr,
                   double* Wrad_vec,
                   double* Wvol_vec,
                   const size_t numParticles,
                   const size_t numDimensions);
*/
#endif //INTERACTION_UTILS_CU_H
