// this experimental version tries to allocate more in the cache size. This is to check how well it will behave for large transmitter arrays.
// since there is excellent locality in cache access, and the cache is prefetched, this may have a minimal impact

// version 5: try to resolve the first line problem
// the problem is that pixels (1,2) (1,3) and so on do not compute correctly. They appear to be shifted 1px to the right - why?
#include "mex.h"
#include "cuda.h"
#include "cuda_runtime.h"
#include "omp.h"
#define pi1 = 3.141592653589793f

static int persistIsInited=0;
static int persistent_deviceCount=0;
/*
 * static int persistent_verbosemode=0;
 * static int persistent_debugvariable1=0;
 */

void beamsim_lambert ( float *tx, unsigned int tx_length, float *out,
        unsigned int n, float d, float r,
        float k, float *xp, float *yp, float *zp)
        // n - no. of pixels //d - density
{
    unsigned int offset=0;
    int ix,iy,itx=0;
    float pressure,distance,kd,pressure_re,pressure_im=0;
    float local_pi=3.141592653589793f;
    float dist2=0;
    float dix,diy,diz,lambert_x,lambert_y,lambert_z=0;
    float xbase,ybase,rho2,rhoi,cosphi,sinphi,cosl,sinl=0;
    float xbase0=-sqrtf((float)(2))+(float)1e-8;
    // di* - delta distances as optimalisation
    #pragma omp parallel for private(iy,ix,pressure_re,pressure_im,xbase,ybase,rho2,rhoi,cosl,cosphi,lambert_x,sinl,lambert_y,sinphi,lambert_z,itx,dix,diy,diz,distance,kd,dist2,pressure,offset) shared(n,tx,xp,yp,zp,out,xbase0,r,d,tx_length,k,local_pi) default(none)
    for (iy=0; iy<n; iy++)
        for (ix=0; ix<n; ix++)
        {
        pressure_re=0;
        pressure_im=0;
        offset=ix+n*iy;
        xbase=(float)ix*d+xbase0;
        ybase=(float)iy*d+xbase0; // it would be an optimisation not to recalculate it for each pixel, this has to stay here be due to future port to CUDA where each pixel has it's own thread
        // this time do exactly the same, but rearrange equations to lower register pressure
        //xxbase=;
        //yybase=; // no need for these
        rho2=xbase*xbase+ybase*ybase; //xxbase and yybase discarded
        if (rho2>2) // optimised from lambert_z<0
        {
            xp[offset]=0;
            yp[offset]=0;
            zp[offset]=0;
            
            out[offset]=(float)0;
            continue;
        }
        rhoi=1/sqrtf(rho2); // rsqrtf = 1/sqrtf but at the cost of a single operation
        cosl=-ybase*rhoi; // rhoi discarded at this point
        cosphi=sqrtf(rho2-rho2*rho2/(float)4);
        lambert_x=r*cosl*cosphi; // cosl discarded at this point
        sinl=xbase*rhoi;
        lambert_y=r*sinl*cosphi; // sinl,cosphi discarded at this point
        
        sinphi=(float)1-rho2/(float)2;
        lambert_z=r*sinphi;
        
        xp[offset]=lambert_x;
        yp[offset]=lambert_y;
        zp[offset]=lambert_z;
        
        for (itx=0; itx<tx_length*6; itx=itx+6)
        {
            dix=(lambert_x-tx[0+itx]);
            diy=(lambert_y-tx[1+itx]);
            diz=(lambert_z-tx[2+itx]);
            distance=sqrtf( dix*dix + diy*diy + diz*diz );
            kd=-k*distance+tx[5+itx];
            dist2=tx[4+itx]/(2*local_pi*distance);
            pressure_re=pressure_re+cosf(kd)*dist2;
            pressure_im=pressure_im+sinf(kd)*dist2;
        }
        // mem write
        pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);
        out[offset]=(float)pressure;

        }
}

__global__ void BeamsimLambertKernel ( float *tx, unsigned int tx_length, float *out,
        unsigned int n, float d, float r,
        float k,float *xp, float *yp, float *zp)
{
    unsigned int offset=0;
    unsigned int ix,iy,itx=0;
    float pressure,distance,kd,pressure_re,pressure_im=0;
    float dist2=0;
    float dix,diy,diz,lambert_x,lambert_y,lambert_z=0;
    float xbase,ybase,rho2,rhoi,cosphi,sinphi,cosl,sinl=0;
    float xbase0=-sqrtf((float)2)+(float)1e-8;
    // calculate ix,iy from thread built-in variables
    ix = blockIdx.x * blockDim.x + threadIdx.x;
    iy = blockIdx.y * blockDim.y + threadIdx.y;
    //ix=0; // debug //
    //C// for (iy=0; iy<ny; iy++)
    //C//  for (ix=0; ix<nx; ix++)
    
    // make sure that this thread won't try to calculate non-existent receiver
    if (iy>n) return;
    
    if (ix>n) return;
    
    // start actual calculation
    pressure_re=0;
    pressure_im=0;
    
    xbase=(float)ix*d+xbase0;
    ybase=(float)iy*d+xbase0; // it would be an optimisation not to recalculate it for each pixel, this has to stay here be due to future port to CUDA where each pixel has it's own thread
    rho2=xbase*xbase+ybase*ybase;
    offset=ix+n*iy;
    if (rho2>(float)2)
    {
        out[offset]=0;
        xp[offset]=0;
        yp[offset]=0;
        zp[offset]=0;
        return;
    }
    rhoi=rsqrtf(rho2);
    cosl=-ybase*rhoi;
    cosphi=sqrtf(rho2-rho2*rho2/(float)4);
    lambert_x=r*cosl*cosphi;
    
    sinl=xbase*rhoi;
    lambert_y=r*sinl*cosphi;
    sinphi=(float)1-rho2/(float)2;
    lambert_z=r*sinphi;
    
    xp[offset]=lambert_x;
    yp[offset]=lambert_y;
    zp[offset]=lambert_z;
    
    for (itx=0; itx<tx_length*6; itx=itx+6) // this hopefully acesses the same memory location for each thread and therefore will be cached
    {
        //    distance=single(sqrt( (ix*dx+x0-tx(1,itx)).^2 + (iy*dy+y0-tx(2,itx)).^2 + (iz*dz+z0-tx(3,itx)).^2 ));
        dix=(lambert_x-tx[0+itx]);
        diy=(lambert_y-tx[1+itx]);
        diz=(lambert_z-tx[2+itx]);
        distance=sqrtf( dix*dix + diy*diy + diz*diz );
        // alternative version:
        //distance=hypotf(hypotf(dix,diy),diz); // turns out that this slows down the code nearly 2x while not providing better accuracy for all realistic parameter values.
        //mexPrintf("distance = %0.2f\n",distance);
        // DirectivityCos=single((iz*dz+z0-tx(3,itx)./distance));
        //C//OBSOLETE// directivitycos=diz/distance;
        // mexPrintf("DirCos = %0.3f\n",directivitycos);
        // kd=(-k*distance+tx(6,itx));
        kd=-k*distance+tx[5+itx];
        //mexPrintf("kd = %0.3f\n",kd);
        // pressure_re=pressure_re+cos(kd)*tx(5,itx)*DirectivityCos/(2*pi*distance);
        //pressure_im=pressure_im+sin(kd)*tx(5,itx)*DirectivityCos/(2*pi*distance);
        
        // oryginal includes directivitycos
        // dist2=tx[4+itx]*directivitycos/(2*pi*distance);
        // *exclude directivitycos, assume ominidirectional until i settle the element directivity
        dist2=tx[4+itx]/(6.283185307179586f*distance); //equals 2*pi
        
//             tmp=tx[4+itx];
//             mexPrintf("tx_amp=%0.2f\n",tmp);
        
        //           mexPrintf("dist2 = %0.3f\n",dist2);
        pressure_re=pressure_re+__cosf(kd)*dist2;
        //          mexPrintf("p_re=%0.3f\n",pressure_re);
        pressure_im=pressure_im+__sinf(kd)*dist2;
        
        // note: __sinf is an simlpified sin function that yields less accurate result. May need to switch to full sin for final product, ok for testing for now
        
        // note 2: function sincosf(...) may be faster in this case - calculates both sin and cos. but since it requires additional accumulators,
        // a detailed test will be required to find out what's faster.
        
    }
    // mem write
    // out(ix+1,iy+1,iz+1)=abs(pressure_re+1i*pressure_im);
    pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);
    // in CUDA, i need to calculate rx array memory offset manually for each thread:
    //offset=ix+nx*iy+(ny*nx)*iz;
    
    // mexPrintf("x = %d; y=%d; z=%d; offset=%d\n",ix,iy,iz,offset);
    out[offset]=(float)pressure; //left for debug
    //debug
    //out[offset]=(float)iy;
    //C// offset=offset++; //this shortcut only works for C
}

static void mexExitFunctionHere()
{
    mexPrintf("beamsim_lambert: exit OK\n");
}

// ----------------- the MEX driver runs on the CPU --------------------
void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    
    //C//float *tx; // 5*n matrix: x,y,z,amplitude,phase
    //C//float *out; // m*
    float k,density,d,r,npts;
    unsigned int n;
    int need_coordinates=0; // return coordinates
    cudaError_t error = cudaSuccess;
    int tmp_ct;
    //C//int totalsize;
    unsigned int tx_length;
    mwSize ndim;
    mwSize dims[2];
    /*  check for proper number of arguments */
    /* NOTE: You do not need an else statement when using mexErrMsgTxt
     * within an if statement, because it will never get to the else
     * statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     * the MEX-file) */
    mexAtExit(mexExitFunctionHere);
    if (persistIsInited==0)
    {
        mexPrintf("beamsim_lambert: Starting beamsim_lambert release 6 kernel. Jerzy Dziewierz, CUE 2008-2013\n");
        // initialize persistent_multidevice_FMCdata to nulls
        cudaGetDeviceCount(&persistent_deviceCount);
        mexPrintf("beamsim_lambert: %d GPUs detected.\n",persistent_deviceCount);
        mexPrintf("beamsim_lambert: This version uses 1 GPU only\n");
        if(persistent_deviceCount == 0){
            mexPrintf("beamsim_lambert: No GPU available, using OpenMP kernel\n");
            int nthreads=omp_get_max_threads();
            mexPrintf("beamsim_lambert: OpenMP: There are %d threads available\n",nthreads);
            if (nthreads<2)
            {
                mexPrintf("beamsim_lambert: Try setting OMP_NUM_THREADS enviroment variable\n");
                mexPrintf("beamsim_lambert: For example, setenv('OMP_NUM_THREADS','2')\n");
            }
            
        }
        persistIsInited=1;
    }
    
    if(nrhs!=4)
        mexErrMsgTxt("4 inputs required.");
    
    if ((nlhs!=1)&&(nlhs!=4))
        mexErrMsgTxt("One or 4 outputs required.");
    if (nlhs==4)
        need_coordinates=1;
    
    
    /* check to make sure the inputs 2-4 input argument are scalar singles */
    for (tmp_ct=1; tmp_ct<4; tmp_ct++)
        if( !mxIsSingle(prhs[tmp_ct]) || mxIsComplex(prhs[tmp_ct]) ||
            mxGetN(prhs[tmp_ct]) * mxGetM(prhs[tmp_ct])!=1 )
        {
        mexErrMsgTxt("Inputs 2-4 must be a scalar of class 'single'.");
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
        mexErrMsgTxt("Tx table is limited to (8*1024*1024) entries for CUDA atm.");
    k=(float)mxGetScalar(prhs[1]);
    r=(float)mxGetScalar(prhs[2]);
    density=(float)mxGetScalar(prhs[3]);
    npts=ceilf(6.283185307179586*r/density);
    d=2*sqrtf((float)2)/npts; // distance between pixels in lambert map
    n=1+(unsigned int)ceilf(2*sqrtf(2)/d); // for some reason this fails to match with pure-C version if i don't add 1 here
    
    ndim=2;
    dims[0]=n; dims[1]=n;
    
    plhs[0] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
    
    mxArray* matlab_xp=mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
    mxArray* matlab_yp=mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
    mxArray* matlab_zp=mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
    if (need_coordinates)
    {
        plhs[1] = matlab_xp;
        plhs[2] = matlab_yp;
        plhs[3] = matlab_zp;
    }
    
    if (persistent_deviceCount>0)
    {
        int tx_arraySize = tx_length * 6 ;
        int tx_memSize = sizeof(float) * tx_arraySize;
        int rx_arraySize = n*n;
        int rx_memSize = sizeof(float) * rx_arraySize;
        
        float *d_tx;
        error=cudaMalloc( &d_tx, tx_memSize );
        if ( error != cudaSuccess )
        {
            mexPrintf("tx_memSize = %d\n",tx_memSize);
            mexPrintf("Error is : %s\n", cudaGetErrorString(cudaGetLastError()));
            mexErrMsgTxt("Memory allocating failure on the GPU.Tx array");
        }
        // plug in the tx array
        
        error=cudaMemcpy( d_tx, (float*) mxGetData(prhs[0]), tx_memSize, cudaMemcpyHostToDevice);
        if (error!=cudaSuccess)
        {
            mexPrintf("Error is : %s\n", cudaGetErrorString(cudaGetLastError()));
            mexErrMsgTxt("Tx Memory copy problem");
        }
        
        
        // now allocate rx array
        
        float *d_rx;
        if ( cudaMalloc( &d_rx, rx_memSize )  != cudaSuccess )
            mexErrMsgTxt("Memory allocating failure on the GPU. Rx array");
        // run kernel
        // assume nz=1 and use 2D blocks, 2D grid
        
        // allocate x,y,z coordinate array
        float *d_xp;
        float *d_yp;
        float *d_zp;
        error=cudaMalloc( &d_xp, rx_memSize );
        error=cudaMalloc( &d_yp, rx_memSize );
        error=cudaMalloc( &d_zp, rx_memSize );
        if ( error != cudaSuccess )
        {
            mexPrintf("rx_memSize = %d\n",rx_memSize);
            mexPrintf("Error is : %s\n", cudaGetErrorString(cudaGetLastError()));
            mexErrMsgTxt("Memory allocating failure on the GPU.Tx array");
        }
        
        dim3 threadsPerBlock(24,16); // threads in a block. 256 threads makes a typical occupancy. reduce if register pressure too high
        dim3 numBlocks((n+threadsPerBlock.x-1)/threadsPerBlock.x,(n+threadsPerBlock.y-1)/threadsPerBlock.y);  // blocks in a grid
        BeamsimLambertKernel<<< numBlocks, threadsPerBlock >>>( d_tx,tx_length,d_rx,n,d,r,k,d_xp,d_yp,d_zp);
        
        /* Get results back from the GPU and free device memory */
        cudaMemcpy( (float*) mxGetData(plhs[0]), d_rx, rx_memSize, cudaMemcpyDeviceToHost);
        if (need_coordinates)
        {
            cudaMemcpy( (float*) mxGetData(plhs[1]), d_xp, rx_memSize, cudaMemcpyDeviceToHost);
            cudaMemcpy( (float*) mxGetData(plhs[2]), d_yp, rx_memSize, cudaMemcpyDeviceToHost);
            cudaMemcpy( (float*) mxGetData(plhs[3]), d_zp, rx_memSize, cudaMemcpyDeviceToHost);
        }
        // finally
        cudaFree( d_tx );
        cudaFree( d_rx );
        cudaFree( d_xp );
        cudaFree( d_yp );
        cudaFree( d_zp );
    }
    else
    { // run the no-cuda version
        float *tx; // 6*n matrix: x,y,z,el_size,amplitude,phase
        float *out; // m*
        float *xp;
        float *yp;
        float *zp;
        tx = (float *)mxGetPr(prhs[0]);
        
        out = (float *)mxGetPr(plhs[0]);
        xp = (float *)mxGetPr(matlab_xp);
        yp = (float *)mxGetPr(matlab_yp);
        zp = (float *)mxGetPr(matlab_zp);
        
        if (need_coordinates)
        {         
        }
        
        beamsim_lambert(tx,tx_length,out,n,d,r,k,xp,yp,zp);
        
        if (need_coordinates==0)
        {         
            // release memory taken by xp,yp,zp, not needed
            mxDestroyArray(matlab_xp);
            mxDestroyArray(matlab_yp);
            mxDestroyArray(matlab_zp);
            
        }
    }
}