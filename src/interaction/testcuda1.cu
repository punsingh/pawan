#include "testcuda1.h"

typedef struct { float4 *pos, *vel; } BodySystem;

void randomizeBodies(float *data, int n) {
    for (int i = 0; i < n; i++) {
        data[i] = 2.0f * (rand() / (float)RAND_MAX) - 1.0f;
    }
}

__global__
void bodyForce(float4 *p, float4 *v, float dt, int n) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < n) {
        float Fx = 0.0f; float Fy = 0.0f; float Fz = 0.0f;

        for (int tile = 0; tile < gridDim.x; tile++) {
            __shared__ float3 spos[BLOCK_SIZE];
            float4 tpos = p[tile * blockDim.x + threadIdx.x];
            spos[threadIdx.x] = make_float3(tpos.x, tpos.y, tpos.z);
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

        v[i].x += dt*Fx; v[i].y += dt*Fy; v[i].z += dt*Fz;
    }
}

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

void testcuda_call(){

    printf("---------------------Entering testcuda_call()----------------------------");
    int num_gpu = 0;  // number of CUDA GPUs
    printf("Launching CUDA computation... \n\n");
    cudaGetDeviceCount(&num_gpu); //get number of gpus available
    if (num_gpu < 1) {
        printf("no CUDA capable GPUs were detected \n");
        return ;
    } else {
        printf("%d CUDA capable GPUs were detected \n", num_gpu);
    }
    int gpuID = num_gpu-1; //the last (non-default) gpu is 'usually' free
    if (cudaSetDevice(gpuID) != cudaSuccess)
        printf("something went wrong setting gpu num %d \n", gpuID);
    cudaDeviceProp cudprop;
    //cudaDeviceProperties_print(&cudprop, gpuID);


    printf("\n \n ");

    int nIters = 1;
    int nBodies;

    for (int iter = 1; iter <= nIters; iter++) {
        nBodies = pow(2,20);
        const float dt = 0.01f; // time step
        printf("nBodies                   = %d \n", nBodies);

        int bytes = 2*nBodies*sizeof(float4);
        std::cout << bytes << std::endl;
        float *buf = (float*)malloc(bytes);
        BodySystem p = { (float4*)buf, ((float4*)buf) + nBodies };

        randomizeBodies(buf, 8*nBodies); // Init pos / vel data

        float *d_buf;
        cudaMalloc(&d_buf, bytes);
        BodySystem d_p = { (float4*)d_buf, ((float4*)d_buf) + nBodies };

        int nBlocks = (nBodies + BLOCK_SIZE - 1) / BLOCK_SIZE;
        printf("No. of BLOCKS               = %d \n No. of THREADS per BLOCK = %d \n", nBlocks, BLOCK_SIZE);

        cudaMemcpy(d_buf, buf, bytes, cudaMemcpyHostToDevice);
        bodyForce<<<nBlocks, BLOCK_SIZE>>>(d_p.pos, d_p.vel, dt, nBodies);
        cudaMemcpy(buf, d_buf, bytes, cudaMemcpyDeviceToHost);
        printf("CUDA compilation finished... \n");
        for (int i = 0 ; i < nBodies; i++) { // integrate position
            p.pos[i].x += p.vel[i].x*dt;
            p.pos[i].y += p.vel[i].y*dt;
            p.pos[i].z += p.vel[i].z*dt;
        }

        free(buf);
        cudaFree(d_buf);
        printf("--------------------- \n");
    }
}
