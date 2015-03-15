// this experimental version tries to allocate more in the cache size. This is to check how well it will behave for large transmitter arrays.
// since there is excellent locality in cache access, and the cache is prefetched, this may have a minimal impact

// version 5: try to resolve the first line problem
// the problem is that pixels (1,2) (1,3) and so on do not compute correctly. They appear to be shifted 1px to the right - why?
#include "mex.h"
#include "cuda.h"
#include "cuda_runtime.h"

//#define pi = 3.141592653589793f
// p2: try to use global constant memory

// -------------- the kernel runs on the GPU ---------------------------

//__device__ __constant__ float tx_const[1536];  // tx table limited to 256 entries - takes 6*256=1536 floats, which is 6144 bytes out of 8192 available
//  __device__ 
        __constant__ float tx_const_cup5[12288]; // allocate 64x64 transmitters
                                   // this is a global variable for GPU so each script has to have a different name, otherwise a collision occurs and "clear all" in matlab is neccesary
// i have noted that this line above actually depends on the driver version. The older driver requires NOT to put in the __device__ in, the newer does.
//cudaMemcpyToSymbol( "dc_ArraySize", &nArr, sizeof(unsigned long) );        

__global__ void timesTwoKernel( float *d_Operand1, float *d_Result1, int arraySize ) {
    int idx = blockDim.x * blockIdx.x + threadIdx.x;
    if (idx > arraySize) return;
    d_Result1[idx] = d_Operand1[idx] * d_Operand1[idx];
} 

__global__ void BeamSimKernel ( float *tx, unsigned int tx_length, float *out,
        float x0, float y0, float z0,
        unsigned int nx, unsigned int ny, unsigned int nz, 
        float dx, float dy, float dz,
        float k)
{
    unsigned int offset=0;
    unsigned int ix,iy,iz,itx=0;
    float pressure,distance,kd,pressure_re,pressure_im=0;
    //C// float directivitycos;
    //float pi=3.141592653589793f;
	float pixel_x,pixel_y,pixel_z,dix,diy,diz,dist2=0;
    float sintmp, costmp;
    //C// float tmp = 0; // not used in CUDA
    
     // di* - delta distances as optimalisation
    
    // iz=0 for now in CUDA, use 2D calculation only
    iz=0;
    
	// C//for (iz=0; iz<nz; iz++)
    // calculate ix,iy from thread built-in variables
    ix = blockIdx.x * blockDim.x + threadIdx.x; 
    // iy=0;
    
    iy = blockIdx.y * blockDim.y + threadIdx.y;
    //ix=0; // debug // 
	//C// for (iy=0; iy<ny; iy++)
    //C//  for (ix=0; ix<nx; ix++)             
    
    // make sure that this thread won't try to calculate non-existent receiver
    if (iy>=ny) return;
    
    if (ix>=nx) return;
    
    // start actual calculation
        pressure_re=0;
        pressure_im=0;

		pixel_x=(float)ix*dx+x0;
		pixel_y=(float)iy*dy+y0; // it would be an optimisation not to recalculate it for each pixel, this has to stay here be due to future port to CUDA where each pixel has it's own thread
		pixel_z=(float)iz*dz+z0;
    
        for (itx=0; itx<tx_length*6; itx=itx+6) // this hopefully acesses the same memory location for each thread and therefore will be cached
         {
            //    distance=single(sqrt( (ix*dx+x0-tx(1,itx)).^2 + (iy*dy+y0-tx(2,itx)).^2 + (iz*dz+z0-tx(3,itx)).^2 ));
            dix=(pixel_x-tx_const_cup5[0+itx]);
            diy=(pixel_y-tx_const_cup5[1+itx]);
            diz=(pixel_z-tx_const_cup5[2+itx]);
            distance=sqrtf( dix*dix + diy*diy + diz*diz ); 
            // alternative version:
            //distance=hypotf(hypotf(dix,diy),diz); // turns out that this slows down the code nearly 2x while not providing better accuracy for all realistic parameter values.
            //mexPrintf("distance = %0.2f\n",distance);
            // DirectivityCos=single((iz*dz+z0-tx(3,itx)./distance));                            
            //C//OBSOLETE// directivitycos=diz/distance;
            // mexPrintf("DirCos = %0.3f\n",directivitycos);
            // kd=(-k*distance+tx(6,itx));
            kd=-k*distance+tx_const_cup5[5+itx];
            //mexPrintf("kd = %0.3f\n",kd);
            // pressure_re=pressure_re+cos(kd)*tx(5,itx)*DirectivityCos/(2*pi*distance);                
             //pressure_im=pressure_im+sin(kd)*tx(5,itx)*DirectivityCos/(2*pi*distance);
            
            // oryginal includes directivitycos
            // dist2=tx[4+itx]*directivitycos/(2*pi*distance);
            // *exclude directivitycos, assume ominidirectional until i settle the element directivity 
            dist2=tx_const_cup5[4+itx]/(6.28318530717958f*distance);
            
//             tmp=tx[4+itx];
//             mexPrintf("tx_amp=%0.2f\n",tmp);
            
 //           mexPrintf("dist2 = %0.3f\n",dist2);
            __sincos(kd,&sintmp,&costmp);
            pressure_re=pressure_re+sintmp*dist2;                    
  //          mexPrintf("p_re=%0.3f\n",pressure_re);
            pressure_im=pressure_im+costmp*dist2; 
            
            // note: __sinf is an simlpified sin function that yields less accurate result. May need to switch to full sin for final product, ok for testing for now
            
            // note 2: function sincosf(...) may be faster in this case - calculates both sin and cos. but since it requires additional accumulators, 
            // a detailed test will be required to find out what's faster.
            
          }        
        // mem write
        // out(ix+1,iy+1,iz+1)=abs(pressure_re+1i*pressure_im); 
        pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);        
        // in CUDA, i need to calculate rx array memory offset manually for each thread:
	    //offset=ix+nx*iy+(ny*nx)*iz;        
        offset=ix+nx*iy;
       // mexPrintf("x = %d; y=%d; z=%d; offset=%d\n",ix,iy,iz,offset);
       out[offset]=(float)pressure; //left for debug
        //debug
       //out[offset]=(float)iy;
       //C// offset=offset++; //this shortcut only works for C                            
}


// ----------------- the MEX driver runs on the CPU --------------------
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
  //C//float *tx; // 5*n matrix: x,y,z,amplitude,phase
  //C//float *out; // m*
  float k,x0,y0,z0,dx,dy,dz;
  unsigned int nx,ny,nz; 
  cudaError_t error = cudaSuccess;  
  int tmp_ct;
  //C//int totalsize;
  unsigned int tx_length;
  mwSize ndim;
  mwSize dims[3];
   /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
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
  if (tx_length>2048)
      mexErrMsgTxt("Tx table is limited to 2048 entries for CUDA atm.");
          
  k=(float)mxGetScalar(prhs[1]);
  x0=(float)mxGetScalar(prhs[2]);
  y0=(float)mxGetScalar(prhs[3]);
  z0=(float)mxGetScalar(prhs[4]);
  nx=(unsigned int)mxGetScalar(prhs[5]);
  ny=(unsigned int)mxGetScalar(prhs[6]);
  nz=(unsigned int)mxGetScalar(prhs[7]);
  if (nz>1)
      mexErrMsgTxt("nz>1 not supported for CUDA for now");
  
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
         mxClassID classid, mxComplexity ComplexFlag */
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
  float *d_tx;
    if ( cudaMalloc( &d_tx, tx_memSize ) != cudaSuccess )
        mexErrMsgTxt("Memory allocating failure on the GPU.Tx array");
  // plug in the tx array
  //OldMethod//
  if (
   cudaMemcpy( d_tx, (float*) mxGetData(prhs[0]), tx_memSize, cudaMemcpyHostToDevice)!=cudaSuccess)
     mexErrMsgTxt("Memory copy problem");
     
  //OldMethod - keep this for now, may not be needed later. Keep for now so that i don't need to change call pattern//
  
  error=cudaMemcpyToSymbol( "tx_const_cup5", (float*) mxGetData(prhs[0]), tx_arraySize*sizeof(float) );
  if(error != cudaSuccess){
		mexPrintf("Error is : %s\n", cudaGetErrorString(cudaGetLastError()));
		mexErrMsgTxt("Problem transferring tx_const table\n");
	}

  
  // now allocate rx array
  int rx_arraySize = nx * ny * nz;
  int rx_memSize = sizeof(float) * rx_arraySize;
  float *d_rx;
  if ( cudaMalloc( &d_rx, rx_memSize )  != cudaSuccess )
        mexErrMsgTxt("Memory allocating failure on the GPU. Rx array");
  // run kernel
  // assume nz=1 and use 2D blocks, 2D grid
  
  dim3 threadsPerBlock(16,24); // threads in a block. 256 threads makes a typical occupancy. reduce if register pressure too high
  // 16,26 is optimal for Fermi. Use 16,24 for earlier cards.
  dim3 numBlocks((nx+threadsPerBlock.x-1)/threadsPerBlock.x,(ny+threadsPerBlock.y-1)/threadsPerBlock.y);  // blocks in a grid
  BeamSimKernel<<< numBlocks, threadsPerBlock >>>( d_tx,tx_length,d_rx,x0,y0,z0,nx,ny,nz,dx,dy,dz,k);
   
  /* Get results back from the GPU and free device memory */
  cudaMemcpy( (float*) mxGetData(plhs[0]), d_rx, rx_memSize, cudaMemcpyDeviceToHost); 
  // finally
  cudaFree( d_tx );
  cudaFree( d_rx );
}