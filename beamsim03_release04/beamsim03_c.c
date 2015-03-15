#include "mex.h"
#include "matrix.h"
#include "math.h"
/*
 * beamsim_proto02.c
 * based on xtimesy.c - example found in API guide
 * direct counterpart to beamsim02_m.m for validation
 * % this version disregards directivity - treats tx as unidirectional points.
*/
/* function out=beamsim02_m(tx,k,x0,y0,z0,nx,ny,nz,dx,dy,dz) */
void beamsim ( float *tx, unsigned int tx_length, float *out,
        float x0, float y0, float z0,
        unsigned int nx, unsigned int ny, unsigned int nz, 
        float dx, float dy, float dz,
        float k)
{
    unsigned int offset=0;
    unsigned int ix,iy,iz,itx=0;
    float pressure,distance,directivitycos,kd,pressure_re,pressure_im=0;
    float pi=3.141592653589793;
	float pixel_x,pixel_y,pixel_z,dix,diy,diz,dist2=0;
    float tmp = 0;
     // di* - delta distances as optimalisation
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
            //    distance=single(sqrt( (ix*dx+x0-tx(1,itx)).^2 + (iy*dy+y0-tx(2,itx)).^2 + (iz*dz+z0-tx(3,itx)).^2 ));
            dix=(pixel_x-tx[0+itx]);
            diy=(pixel_y-tx[1+itx]);
            diz=(pixel_z-tx[2+itx]);
            distance=sqrtf( dix*dix + diy*diy + diz*diz ); 
            //mexPrintf("distance = %0.2f\n",distance);
            // DirectivityCos=single((iz*dz+z0-tx(3,itx)./distance));                
            directivitycos=diz/distance;
            // mexPrintf("DirCos = %0.3f\n",directivitycos);
            // kd=(-k*distance+tx(6,itx));
            kd=-k*distance+tx[5+itx];
            //mexPrintf("kd = %0.3f\n",kd);
            // pressure_re=pressure_re+cos(kd)*tx(5,itx)*DirectivityCos/(2*pi*distance);                
             //pressure_im=pressure_im+sin(kd)*tx(5,itx)*DirectivityCos/(2*pi*distance);
            
            // oryginal includes directivitycos
            // dist2=tx[4+itx]*directivitycos/(2*pi*distance);
            // *exclude directivitycos, assume ominidirectional until i settle the element directivity 
            dist2=tx[4+itx]/(2*pi*distance);
            
//             tmp=tx[4+itx];
//             mexPrintf("tx_amp=%0.2f\n",tmp);
            
 //           mexPrintf("dist2 = %0.3f\n",dist2);
            pressure_re=pressure_re+cosf(kd)*dist2;                    
  //          mexPrintf("p_re=%0.3f\n",pressure_re);
            pressure_im=pressure_im+sinf(kd)*dist2; 
            // skip the imaginary part
          }        
        // mem write
        // out(ix+1,iy+1,iz+1)=abs(pressure_re+1i*pressure_im); 
        pressure=sqrtf(pressure_re*pressure_re+pressure_im*pressure_im);        
	    //offset=ix+nx*iy+(ny*nx)*iz;        
       // mexPrintf("x = %d; y=%d; z=%d; offset=%d\n",ix,iy,iz,offset);
        out[offset]=(float)pressure;
        offset=offset++;
       }
          
    
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    
  float *tx; // 6*n matrix: x,y,z,el_size,amplitude,phase
  float *out; // m*
  float k,x0,y0,z0,dx,dy,dz;
  unsigned int nx,ny,nz; 
    
  int tmp_ct;
  int totalsize;
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
  
  k=(float)mxGetScalar(prhs[1]);
  x0=(float)mxGetScalar(prhs[2]);
  y0=(float)mxGetScalar(prhs[3]);
  z0=(float)mxGetScalar(prhs[4]);
  nx=(unsigned int)mxGetScalar(prhs[5]);
  ny=(unsigned int)mxGetScalar(prhs[6]);
  nz=(unsigned int)mxGetScalar(prhs[7]);
  if (nx>4096)
     mexErrMsgTxt("nx must be a sensible value 0<nx<4096");     
  if (ny>4096)
     mexErrMsgTxt("ny must be a sensible value 0<ny<4096");     
  if (nz>4096)
     mexErrMsgTxt("nz must be a sensible value 0<nz<4096");     
  
  dx=(float)mxGetScalar(prhs[8]);  
  dy=(float)mxGetScalar(prhs[9]);
  dz=(float)mxGetScalar(prhs[10]);

  totalsize = nx*ny*nz;  

  ndim=3;

  dims[0]=nx; dims[1]=ny; dims[2]=nz;
  
  plhs[0] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
          /* mwSize ndim, const mwSize *dims, 
         mxClassID classid, mxComplexity ComplexFlag */
  /*  create a C pointer to a copy of the output matrix */
  out = (float *)mxGetPr(plhs[0]);
  tx = (float *)mxGetPr(prhs[0]);
  //out[0]=(float)1; out[1]=(float)2; 
  beamsim(tx,tx_length,out,x0,y0,z0,nx,ny,nz,dx,dy,dz,k);
  
}
