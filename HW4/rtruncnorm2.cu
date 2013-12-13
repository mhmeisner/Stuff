#include <stdio.h>
#include <stdlib.h>

#include <cuda.h>
#include <curand_kernel.h>
#include <math_constants.h>

extern "C"
{

__global__ void 
rtruncnorm_kernel(
    float *x, int n,
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
        curandState_t rng;
        curand_init(rng_a+rng_b*idx,rng_c,0,&rng);

        // Draw sample
        int ntries = 0;
        int accepted = 0;
        float ran;
        while(!accepted and ntries < maxtries){
            ran = mu[idx]+sigma[idx]*curand_normal(&rng);
            ntries += 1;
            if(ran >= lo[idx] and ran <= hi[idx]){
                accepted = 1;
            }
        }
        // Use Robert Method if that didn't work
        if(!accepted){
            // my code is set up to sample only (a,infty), so if it's a (-infty,b) sample we want, then we sample from (-b,infty) and reverse the sign after
            int rev_sign = 0;
            float lower;
            if(isfinite(hi[idx])){
                lower = lo[idx]-mu[idx];
            }else{
                lower = mu[idx]-hi[idx];
                rev_sign = 1;
            }
            float alpha = (lower+sqrtf(lower*lower+4))/2;
            float z;
            int ntries = 0;
            // I may well have done something wrong...but for some datasets, this while loop never ended if I didn't set a max # of tries.
            while(!accepted and ntries < 10000L){
                ntries += 1;
                float psi;
                // sample uniform, then use inverse cdf to get sample from exponential distribution: 
                z = lower-logf(curand_uniform(&rng))/alpha;
                if(lower<alpha){
                    psi = expf(-powf((alpha-z),2)/2);
                }else{
                    psi = expf(-powf((alpha-z),2)/2)*expf(powf((lower-alpha),2)/2);
                }
                float u = curand_uniform(&rng);
                if(u<psi){
                    accepted = 1;
                }
            }
            if(rev_sign){
                ran = mu[idx]-z;
            }else{
                ran = mu[idx]+z;
            }
            // If the Robert method failed to accept in 10000 tries:
            if(!accepted){
                ran = CUDART_NAN_F;
            }
        }
       
        x[idx] = ran;
    }

    return;
}

} // END extern "C"

