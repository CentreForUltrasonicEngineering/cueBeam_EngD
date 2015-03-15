# cueBeam: acoustic beam simulation using Huygens's principle.

cueBeam is a Matlab/CUDA code for simulation of acoustic field pressure distribution. It calculates pressure in frequency domain (meaning, assumes that the radiators continuously radiate) using the Huygen's principle.

This ultra-simplified propagation model enables obtaining quick estimates of pressure field shape, beam width, and side lobe amplitude, suitable for use in NDT research. It is admittedly less accurate than other published methods, but the advantage is in it's speed.

Radiators are described in terms of spatial XYZ location, amplitude and phase; despite them being point-like, they physically represent perfectly baffled pistons.

Propagation medium is described in terms of wavenumber k. Radiating sources are treated as point like. For each field probing point, value is computed by summing the contributions from each transmitter, shifted by appropriate phase delay and reduced by respective distance between the transmitter and probed point:


    for each pixel
      	pressure=complex zero;
    	for each radiator
    		distance=distance(radiator,pixel);
    		phase_shift=wavenumber*distance+radiator_phaseshift;
    		amplitdue_decayed=radiator_amplitude/distance;
    		pressure=pressure+...
    			complex(amplitdue_decayed,phase_shift);
    	end for each radiator
    	out(pixel)=abs(pressure);  % store absolute value
    end for each pixel

## Preparation

The following inputs are required:

1. Locations of the centre of radiating elements -x,y,z.
2. 'steering vector' - a complex number describing amplitude and phase of radiation for each element. The size of the SVect must be n*1
3. Wavenumber k in the medium (single frequency only)
4. Description of the location of the field probing points. There are two versions for this: planar and spherical. The planar is called XZ (planar sampling space along XZ plane) and the spherical is called Lambert, for it uses equi-areal azimuthal sampling of hemisphere.

 These inputs should be packed into a n*6 matrix

    tx=[elemX elemY elemZ zeros() abs(SVect(:)) angle(SVect(:))];
    

## Planar sampling version

For XZ sampling, the sampling points (where pressure is calculated) is always a regular grid described by:  (x0,y0,z0) - location of a corner of the grid; (dx,dy,dz) - distance between points, and (nx,ny,nz) : number of points in each direction.  ny=1 always in this implementation.

The way to call the calculation function is:

	img_xz = cueBeam_xz(tx',k,x0,y0,z0,nx,ny,nz,dx,dy,dz);

Note that all inputs must be of class single. (default for Matlab is double, so conversion is needed). This is both for speed and compatibility with early CUDA cards. No ill conditioned math is involved so single precision numbers deliver approximately 80dB of precision.

![space](cueBeam_XZ_space.png)
 
> Figure 1. Scene setting for cueBEAM. Green: element locations(actually, points); red cross: "probe points" where field is calculated. Note that this is an old version of the figure and the default probing points are XZ only.(**todo: create updated image showing correct XZ plane!**)

## Hemisphere sampling version

For beamsim Lambert the calling convention is simpler:

	[img_lambert lambert_x lambert_y lambert_z]=cueBeam.cueBeam_lambert(tx', k, r, density); 

This program automatically generates a mesh of points that are distributed over a hemisphere with radius R, and distance between real points "density", in such way, that the true area covered by a given point is equal for all points. The points are then mapped to a rectangle. The transformation rules are (after http://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html )

R is the radius of the sphere, and density is an approximate linear distance between points on the sphere. The resulting resolution of the field image can be calculated by:

    number_of_points=ceil(2*pi*r/density);
    d=2*sqrt(2)/number_of_points;
    n=ceil(2*sqrt(2)/d); % note, the double scaling is needed to obtain correct rounding.
    img_lambert=zeros(n,n,'single'); 
     

Inverse transformation rules are used to calculate standard parallel   and central longitude on the sphere:   
 
![Lambert 1](cueBeam_Lambert_space_1.png)
![Lambert 2](cueBeam_Lambert_space_2.png) 
>Figure 2. Lambert azimuthal equiareal map field probing point distribution.


When compared to regular orthogonal grid, the advantage of this approach is that true power of the sidelobes can be easily integrated and compared with the power of the main lobe. This reduces error of estimating the sidelobe level for classic beamforming, and therefore is more representative of the image contrast from the operator's point of view. 

## Example

An example on how to define parameters for use with cueBeam are provided in the MATLAB file:

cueBeam_linear_array_basic


