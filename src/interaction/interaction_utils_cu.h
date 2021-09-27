//
// Created by Sumeet Kumar on 26/08/21.
//

#ifndef INTERACTION_UTILS_CU_H
#define INTERACTION_UTILS_CU_H


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
