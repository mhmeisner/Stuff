#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{

__global__ void 
rtruncnorm_kernel(float *x, int n, 
                  float *mu, float *sigma, 
                  float *lo, float *hi,
                  int maxtries, int rng_a, 
                  int rng_b, int rng_c)
{
    // Usual block/thread indexing...
    int myblock = blockIdx.x + blockIdx.y * gridDim.x;
    int blocksize = blockDim.x * blockDim.y * blockDim.z;
    int subthread = threadIdx.z*(blockDim.x * blockDim.y) + threadIdx.y*blockDim.x + threadIdx.x;
    int idx = myblock * blocksize + subthread;

    if(idx < n){
        // Setup the RNG:
        curandState rng;
        curand_init(rng_a+rng_b*idx,rng_c,0,&rng);

        // Draw sample
        int ntries = 0;
        int accepted = 0;
        while(!accepted and ntries < maxtries){
            float ran = curand_normal(&rng);
            ran = mu[idx]+ran*sigma[idx];
            ntries += 1;
            if(ran >= lo[idx] and ran <= hi[idx]){
                accepted = 1;
            }
        }

        // Store sample:
        x[idx] = ran;
    }

    return;
}

} // END extern "C"

