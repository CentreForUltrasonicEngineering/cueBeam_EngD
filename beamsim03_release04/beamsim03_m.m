% beamsim03_m 
% Matlab m-script version of the beam shape calculation code
% trivial version designed to be convertible to C and CUDA
%
% Calculates beam shape at single frequency.
% 
% use:
% tx=single([tx_x(:) tx_y(:) tx_z(:) tx_size abs(SVect(:)) angle(SVect(:))])'; 
%
% where 
% -> tx_x,tx_y,tx_z - locations of the transmitters
%
% -> tx_size - size of the transmitter, assuming it's circular. not
% implemented in this version, sits here for compatibility with future
% versions
% 
% -> SVect - Steering vector - a complex number for each transmitter -
% amplitude and phase of the tx signal at given frequency
%
% tx is a (6,:) matrix because of the way matlab stores the data
%
% and then 
% out=beamsim03_m(tx,k,x0,y0,z0,nx,ny,nz,dx,dy,dz)
%
% -> tx - description of the transmitters
% -> k - wavenumber
% -> x0,y0,z0 - point in space where image starts
% -> nx,ny,nz - number of pixels in each dimension
% -> dx,dy,dx - spacing between pixels
%
% Jerzy Dziewierz, CUE 2010
%

function out=beamsim03_m(tx,k,x0,y0,z0,nx,ny,nz,dx,dy,dz)
% this version disregards directivity - treats tx as unidirectional points.

% reserve memory

out=zeros(nx,ny,nz,'single');
%k = 2 * pi / wavelength;
for iz=0:(nz-1)
    for iy=0:(ny-1)
         for ix=0:(nx-1)
            pressure_re=single(0);
            pressure_im=single(0);
            for itx=1:size(tx,2)                
                distance=single(sqrt( (single(ix)*dx+x0-tx(1,itx)).^2 + (single(iy)*dy+y0-tx(2,itx)).^2 + (single(iz)*dz+z0-tx(3,itx)).^2 ));
                kd=(-k*distance+tx(6,itx));
                pressure_re=pressure_re+cos(kd)*tx(5,itx)/(2*pi*distance);
                pressure_im=pressure_im+sin(kd)*tx(5,itx)/(2*pi*distance);
            end
            % write out result for that pixel
            out(ix+1,iy+1,iz+1)=abs(pressure_re+1i*pressure_im); 
        end
    end
end
%% helper
% tx=single([tx_x(:) tx_y(:) tx_z(:) tx_size abs(SVect(:)) angle(SVect(:))])'; 
% tx is a (6,:) matrix because of the way matlab stores the data