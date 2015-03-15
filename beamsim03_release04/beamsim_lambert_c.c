#include "mex.h"
#include "matrix.h"
#include "math.h"
#define pi 3.14159265358979323846f
/*
 * beamsim_0303o2_c.c
 * based beamsim03
 * direct counterpart to beamsim0303o2.m for validation
 * inputs: [out]=beamsim0303o2_m(tx,k,r,d)
 * % this version disregards directivity - treats tx as unidirectional points.
*/
/* function out=beamsim02_m(tx,k,x0,y0,z0,nx,ny,nz,dx,dy,dz) */
void beamsim_lambert ( float *tx, unsigned int tx_length, float *out,
        unsigned int n, float d, float r,
        float k)
        // n - no. of pixels //d - density
{ 
    unsigned int offset=0;
    unsigned int ix,iy,itx=0;
    float pressure,distance,kd,pressure_re,pressure_im=0;
//    float pi=3.141592653589793f;
	float dist2=0;
    float dix,diy,diz,lambert_x,lambert_y,lambert_z=0;
    float xbase,ybase,rho2,rhoi,cosphi,sinphi,cosl,sinl=0;
    float tmp = 0;
    float xbase0=-sqrtf((float)(2))+(float)1e-8;
     // di* - delta distances as optimalisation
	 for (iy=0; iy<n; iy++)
      for (ix=0; ix<n; ix++)             
      {
        pressure_re=0;
        pressure_im=0;

		xbase=(float)ix*d+xbase0;
		ybase=(float)iy*d+xbase0; // it would be an optimisation not to recalculate it for each pixel, this has to stay here be due to future port to CUDA where each pixel has it's own thread		
          // this time do exactly the same, but rearrange equations to lower register pressure
          //xxbase=;
          //yybase=; // no need for these
          rho2=xbase*xbase+ybase*ybase; //xxbase and yybase discarded                             
          if (rho2>2) // optimised from lambert_z<0
          {
              out[offset]=(float)0;
              offset=offset++;
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
          
          
          for (itx=0; itx<tx_length*6; itx=itx+6)
         {
                     
            //    distance=single(sqrt( (ix*dx+x0-tx(1,itx)).^2 + (iy*dy+y0-tx(2,itx)).^2 + (iz*dz+z0-tx(3,itx)).^2 ));
            dix=(lambert_x-tx[0+itx]);
            diy=(lambert_y-tx[1+itx]);
            diz=(lambert_z-tx[2+itx]);
            distance=sqrtf( dix*dix + diy*diy + diz*diz );             
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
  float k,density,d,r;
    
  int tmp_ct;
  int totalsize;
  unsigned int n;
  float npts;
  unsigned int tx_length;
  
  mwSize ndim;
  mwSize dims[2];
   /*  check for proper number of arguments */
  /* NOTE: You do not need an else statement when using mexErrMsgTxt
     within an if statement, because it will never get to the else
     statement if mexErrMsgTxt is executed. (mexErrMsgTxt breaks you out of
     the MEX-file) */
  if(nrhs!=4) 
    mexErrMsgTxt("4 inputs required.");
  if(nlhs!=1) 
    mexErrMsgTxt("One output required.");
  
  
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
  
  k=(float)mxGetScalar(prhs[1]);
  r=(float)mxGetScalar(prhs[2]);
  density=(float)mxGetScalar(prhs[3]);
  
  npts=(float)ceil(2*pi*r/density);
  
  d=2*sqrtf((float)2)/npts; // distance between pixels in lambert map
  
  n=(unsigned int)ceil(2*sqrtf(2)/d);
  totalsize = n*n;

  ndim=2;

  dims[0]=n; dims[1]=n;
  
  plhs[0] = mxCreateNumericArray(ndim,dims,mxSINGLE_CLASS,mxREAL);
          /* mwSize ndim, const mwSize *dims, 
         mxClassID classid, mxComplexity ComplexFlag */
  /*  create a C pointer to a copy of the output matrix */
  out = (float *)mxGetPr(plhs[0]);
  tx = (float *)mxGetPr(prhs[0]);
  //out[0]=(float)1; out[1]=(float)2; 
  
  beamsim_lambert(tx,tx_length,out,n,d,r,k);
  
}
