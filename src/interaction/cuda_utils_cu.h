//
// Created by Sumeet Kumar on 12/09/21.
//

#ifndef CUDA_UTILS_CU_H
#define CUDA_UTILS_CU_H

void cudaDeviceProperties_print(cudaDeviceProp *cudprop, const int &gpuID){
    cudaGetDeviceProperties(cudprop, gpuID);
    printf("-----------device properties------------- \n");
    printf("| \t cudaDeviceProp.name                : \t %s \n", cudprop->name);
    printf("| \t cudaDeviceProp.maxThreadsPerBlock  : \t %d \n", cudprop->maxThreadsPerBlock);
    printf("| cudaDeviceProp.maxThreadsPerMultiProcessor: \t %d \n", cudprop->maxThreadsPerMultiProcessor);
    // printf("| cudaDeviceProp.maxBlocksPerMultiProcessor : \t %d \n", cudprop->maxBlocksPerMultiProcessor);
    printf("| \t cudaDeviceProp.maxThreadsDim       : \t %d,%d,%d \n", cudprop->maxThreadsDim[0],
           cudprop->maxThreadsDim[1],
           cudprop->maxThreadsDim[2]);
    printf("| \t cudaDeviceProp.maxGridSize         : \t %d,%d,%d \n", cudprop->maxGridSize[0],
           cudprop->maxGridSize[1],
           cudprop->maxGridSize[2]);
    printf("| \t cudaDeviceProp.warpSize            : \t %d \n", cudprop->warpSize);
    printf("----------------------------------------- \n");

}


#endif //CUDA_UTILS_CU_H
