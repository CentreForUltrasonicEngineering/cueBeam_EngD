#include "mex.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "omp.h"

// #define FORCE_NOCUDA

#define pi2 2.0f*3.141592653589793f
#define pi1 3.141592653589793f
static int persistIsInited=0;
static int persistent_deviceCount=0;
/*static int persistent_verbosemode=0;
static int persistent_debugvariable1=0;*/

__host__ void beamsim ( float *tx, unsigned int tx_length, float *out,
        float x0, float y0, float z0,
        unsigned int nx, unsigned int ny, unsigned int nz,
        float dx, float dy, float dz,
        float k)
{
    unsigned int offset=0;
    int ix,iy,iz,itx=0;
    float pressure,distance,kd,pressure_re,pressure_im=0;
    float pixel_x,pixel_y,pixel_z,dix,diy,diz,dist2=0;
    // di* - delta distances as optimalisation
    #pragma omp parallel for private(iz,iy,ix,pressure_re,pressure_im,pixel_x,pixel_y,pixel_z,itx,dix,diy,diz,distance,kd,dist2,pressure,offset) shared(ny,nx,tx)
    for (iz=0; iz<nz; iz++)        
        for (iy=0; iy<ny; iy++)
            for (ix=0; ix<nx; ix++)
            {
        pressure_re=0;
        pressure_im=0;
        
        pixel_x=(float)ix*dx+x0;
        pixel_y=(float)iy*dy+y0; // it would be an optimisation not to recalculate it for each pixel, this has to stay here be due to future port to CUDA where each pixel has it's own thread
        pixel_z=(float)iz*dz+z0;
        
        for (itx=0; itx<tx_length*6; itx=itx+6)
        {
            dix=(pixel_x-tx[0+itx]);
            diy=(pixel_y-tx[1+itx]);
            diz=(pixel_z-tx[2+itx]);
            distance=sqrtf( dix*dix + diy*diy + diz*diz );
            kd=-k*distance+tx[5+itx];
            dist2=tx[4+itx]/(2*pi1*distance);
            pressure_re=pressure_re+cosf(kd)*dist2;
            pressure_im=pressure_im+sinf(kd)*dist2;
            
        }
        // mem write
        pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);
        offset=ix+nx*iz;
        out[offset]=(float)pressure;
            }
}

__global__ void BeamSimKernel ( const float *tx, unsigned int tx_length, float *out,
        float x0, float y0, float z0,
        unsigned int nx, unsigned int ny, unsigned int nz,
        float dx, float dy, float dz,
        float k)
{
    unsigned int offset=0;    
    unsigned int ix,iy,iz,itx=0;    
    float pressure,distance,kd,pressure_re,pressure_im=0;    
    float pixel_x,pixel_y,pixel_z,dix,diy,diz,dist2=0;    
    
    iz=blockIdx.y * blockDim.y + threadIdx.y;; // z is used on y axis of the device    
    ix = blockIdx.x * blockDim.x + threadIdx.x;    
    iy = 0;    
    
    // make sure that this thread won't try to calculate non-existent receiver
    if (iy>=ny) return;    
    if (ix>=nx) return;    
    if (iz>=nz) return;
    
    // start actual calculation
    pressure_re=0;
    pressure_im=0;
    
    pixel_x=(float)ix*dx+x0;
    pixel_y=(float)iy*dy+y0;
    pixel_z=(float)iz*dz+z0;
    
    for (itx=0; itx<tx_length*6; itx=itx+6) // this hopefully acesses the same memory location for each thread and therefore will be cached
    {        
        dix=(pixel_x-tx[0+itx]);
        diy=(pixel_y-tx[1+itx]);
        diz=(pixel_z-tx[2+itx]);
        distance=sqrtf( dix*dix + diy*diy + diz*diz );        
        kd=-k*distance+tx[5+itx];        
        dist2=tx[4+itx]/(pi2*distance);
        pressure_re=pressure_re+__cosf(kd)*dist2;
        pressure_im=pressure_im+__sinf(kd)*dist2;
    }    
    // mem write
    pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);    
    offset=ix+nx*iz;    
    out[offset]=(float)pressure; //left for debug
}

static void mexExitFunctionHere()
{
    mexPrintf("beamsim_xz: exit OK\n");
}

// ----------------- the MEX driver runs on the CPU --------------------
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    float k,x0,y0,z0,dx,dy,dz;
    unsigned int nx,ny,nz;
    cudaError_t error = cudaSuccess;
    int tmp_ct;
    unsigned int tx_length;
    mwSize ndim;
    mwSize dims[3];
    /*  check for proper number of arguments */
    mexAtExit(mexExitFunctionHere);
    
    if (persistIsInited==0)
    {
        mexPrintf("beamsim_xz: Starting beamsim_xz release 6 kernel. Jerzy Dziewierz, CUE 2008-2013\n");
        // initialize persistent_multidevice_FMCdata to nulls
        cudaGetDeviceCount(&persistent_deviceCount);
        mexPrintf("beamsim_xz: %d GPUs detected.\n",persistent_deviceCount);
        mexPrintf("beamsim_xz: This version uses 1 GPU only\n");
        if(persistent_deviceCount == 0){
            mexPrintf("beamsim_xz: No GPU available, using OpenMP kernel\n");
            int nthreads=-1;
            #pragma omp parallel
            #pragma omp master
            {
                nthreads = omp_get_num_threads();
                mexPrintf("beamsim_xz: OpenMP: There are %d threads available\n",nthreads);
            } 
            
            if (nthreads<2)
            {
                mexPrintf("beamsim_xz: Try setting OMP_NUM_THREADS enviroment variable\n");
                mexPrintf("beamsim_xz: For example, setenv('OMP_NUM_THREADS','2')\n");
            }
            
        }
        persistIsInited=1;
        #ifdef FORCE_NOCUDA                 
        mexPrintf("Note: Forcing NO CUDA\n");
        persistent_deviceCount=0;
        #endif
        
    }
    
    if(nrhs!=11)
        mexErrMsgTxt("11 inputs required.");
    
    if(nlhs!=1)
        mexErrMsgTxt("One output required.");
    
    
    /* check to make sure the inputs 2-5 and 9-11 input argument are scalar singles */
    for (tmp_ct=1; tmp_ct<5; tmp_ct++)
        if( !mxIsSingle(prhs[tmp_ct]) || mxIsComplex(prhs[tmp_ct]) ||
            mxGetN(prhs[tmp_ct]) * mxGetM(prhs[tmp_ct])!=1 ) {
        mexErrMsgTxt("Inputs 2-5 must be a scalar of class 'single'.");
        }
    for (tmp_ct=8; tmp_ct<11; tmp_ct++)
        if( !mxIsSingle(prhs[tmp_ct]) || mxIsComplex(prhs[tmp_ct]) ||
            mxGetN(prhs[tmp_ct]) * mxGetM(prhs[tmp_ct])!=1 ) {
        mexErrMsgTxt("Inputs 9-11 must be a scalar of class 'single'.");
        }
    /* check to make sure that inputs 10,11,12 are scalar integers */
    for (tmp_ct=5; tmp_ct<8; tmp_ct++)
        if( !mxIsUint32(prhs[tmp_ct]) || mxIsComplex(prhs[tmp_ct]) ||
            mxGetN(prhs[tmp_ct]) * mxGetM(prhs[tmp_ct])!=1 ) {
        mexErrMsgTxt("Inputs 6-8 must be a scalar of class 'uint32'.");
        }
    /* check to make sure that input 1 is a n*6 single matrix */
    if (!mxIsSingle(prhs[0]) || mxIsComplex(prhs[0]) ||
            mxGetM(prhs[0])!=6)
    {
        mexErrMsgTxt("Input 1 must be 6*n matrix of class 'single'.");
    }
    
    // get variables from input
    tx_length=(unsigned int)mxGetN(prhs[0]);
    if (tx_length>(8*1024*1024))
        mexErrMsgTxt("Tx table is limited to 8*1024*1024 entries for CUDA atm.");
    
    k=(float)mxGetScalar(prhs[1]);
    x0=(float)mxGetScalar(prhs[2]);
    y0=(float)mxGetScalar(prhs[3]);
    z0=(float)mxGetScalar(prhs[4]);
    nx=(unsigned int)mxGetScalar(prhs[5]);
    ny=(unsigned int)mxGetScalar(prhs[6]);
    nz=(unsigned int)mxGetScalar(prhs[7]);
    
    if (ny>1)
        mexErrMsgTxt("This version only works along XZ axes. Try xy or yz version");
    
    if (nx>4096)
        mexErrMsgTxt("nx must be a sensible value 0<nx<4096");
    if (ny>4096)
        mexErrMsgTxt("ny must be a sensible value 0<ny<4096");
    if (nz>4096)
        mexErrMsgTxt("nz must be a sensible value 0<nz<4096");
    
    dx=(float)mxGetScalar(prhs[8]);
    dy=(float)mxGetScalar(prhs[9]);
    dz=(float)mxGetScalar(prhs[10]);
    
    //C//totalsize = nx*ny*nz;
    
    ndim=3;
    
    dims[0]=nx; dims[1]=ny; dims[2]=nz;
    
    plhs[0] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
    
    /* mwSize ndim, const mwSize *dims,
     * mxClassID classid, mxComplexity ComplexFlag */
    /*  create a C pointer to a copy of the output matrix */
    //C//out = (float *)mxGetPr(plhs[0]);
    //C//tx = (float *)mxGetPr(prhs[0]);
    
    //out[0]=(float)1; out[1]=(float)2;
    // oryginal beamsim
    //beamsim(tx,tx_length,out,x0,y0,z0,nx,ny,nz,dx,dy,dz,k);
    
    // it is now time to allocate memory in gpu
    
    // first goes the tx array
    // this will go into constant memory later but resides in normal mem for now
    int tx_arraySize = tx_length * 6 ;
    int tx_memSize = sizeof(float) * tx_arraySize;
    int rx_arraySize = nx * ny * nz; // nz y is always 1 here
    int rx_memSize = sizeof(float) * rx_arraySize;
    
    // version for CUDA enabled
    if (persistent_deviceCount>0)
    {
        float *d_tx; // device tx table
        error=cudaMalloc( &d_tx, tx_memSize );
        if (  error!= cudaSuccess )
        {
            mexPrintf("Error is : %s\n", cudaGetErrorString(cudaGetLastError()));
            mexErrMsgTxt("Memory allocating failure on the GPU.Tx array");
        }
        // plug in the tx array
        error=cudaMemcpy( d_tx, (float*) mxGetData(prhs[0]), tx_memSize, cudaMemcpyHostToDevice);
        if ( error!=cudaSuccess )
        {
            mexPrintf("Error is : %s\n", cudaGetErrorString(cudaGetLastError()));
            mexErrMsgTxt("Memory copy problem");
        }
        
        
        float *d_rx;
        error=cudaMalloc( &d_rx, rx_memSize );
        if ( error  != cudaSuccess )
        {
            mexPrintf("Error is : %s\n", cudaGetErrorString(cudaGetLastError()));
            mexErrMsgTxt("Memory allocating failure on the GPU. Rx array");
        }
        // run kernel
        
        // assume ny=1 and use 2D blocks, 2D grid. z is used at y axis
        
        dim3 threadsPerBlock(16,26); // threads in a block. 256 threads makes a typical occupancy. reduce if register pressure too high
        // 16,26 is optimised for fermi. Use 16/24 for earlier cards.
        
        dim3 numBlocks((nx+threadsPerBlock.x-1)/threadsPerBlock.x,(nz+threadsPerBlock.y-1)/threadsPerBlock.y);  // blocks in a grid
        
        BeamSimKernel<<< numBlocks, threadsPerBlock >>>( d_tx,tx_length,d_rx,x0,y0,z0,nx,ny,nz,dx,dy,dz,k);
        
        /* Get results back from the GPU and free device memory */
        cudaMemcpy( (float*) mxGetData(plhs[0]), d_rx, rx_memSize, cudaMemcpyDeviceToHost);
        // finally
        cudaFree( d_tx );
        cudaFree( d_rx );
    } // if CUDA available
    else
    { 
        float *out;
        float *tx;
// run no-cuda version       
       out = (float *)mxGetPr(plhs[0]);
       tx = (float *)mxGetPr(prhs[0]);
    beamsim(tx,tx_length,out,x0,y0,z0,nx,ny,nz,dx,dy,dz,k);
    
    }
    
}