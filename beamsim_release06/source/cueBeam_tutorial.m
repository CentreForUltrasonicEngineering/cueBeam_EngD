% beamsim_tutorial
% this script demonstrates how to use the beamsim set.

% set up enviroment parameters. All units in SI
frequency=1e6;
wave_velocity=5600; % steel
wavelength=wave_velocity/frequency;
wavenumber=single(2*pi/wavelength);

% set XYZ of focal point
alpha_rotation=-pi/12;
focal_distance=50e-3;
focal_x=focal_distance*sin(alpha_rotation); focal_y=0; focal_z=focal_distance*cos(alpha_rotation);

% set up probe element locations.
% assume 2D phased array probe with 8*8 elements
array_size=50e-3;
elements_in_row=8;

% generate array element locations based on array_size and elements_in_row
x_row=linspace(-array_size/2,array_size/2,elements_in_row);
y_row=x_row;
[element_x element_y]=meshgrid(x_row,y_row);
element_x=element_x(:); element_y=element_y(:); % convert to linear table
element_z=zeros(size(element_x));
element_amplitude=1+zeros(size(element_x)); % driving amplitude

% calculate phase for each element in such way, that the beam is focussed
% in the focal point

distances = sqrt((element_x-focal_x).^2+(element_y-focal_y).^2+(element_z-focal_z).^2); 
SVect=element_amplitude.*exp(1i*2*pi*distances/wavelength); % steering vector. need to update that in the callback too!

% assemble tx vector
tx=single([element_x element_y element_z zeros(size(element_x)) abs(SVect(:)) angle(SVect(:))]);

% define image parameters
lambert_radius=single(50e-3);
lambert_map_density=single(0.5e-3);

% get lambert map image
tic
img_lambert=cueBeam_lambert(tx',wavenumber,lambert_radius,lambert_map_density);           
tout=toc; raycount=prod(size(img_lambert))*size(tx,1); rayspeed=raycount/tout; fprintf('%0.1f Mrays/s\n',rayspeed/1e6);
% convert image to decibel scale
img_db=20*log10(img_lambert./max(img_lambert(:)));
lambert_coords=[-90 90];
% display lambert map image

figure(1);
imagesc(lambert_coords,lambert_coords,img_db,[-20 0]); colorbar; axis image;
title('Lambert map : dB range');
xlabel('angle[deg]'); ylabel('angle[deg]');

%%
% get XZ cross-section image
% define XZ image parameters
resolution=single(0.05e-3);
x0=single(-20e-3); x1=single(20e-3); dx=resolution; nx=uint32(ceil((x1-x0)/dx));
y0=single(0); y1=single(0); dy=resolution; ny=uint32(ceil((y1-y0)/dy)+1);
z0=single(10e-3); z1=single(60e-3); dz=resolution; nz=uint32(ceil((z1-z0)/dz));
tic;
img_xz=squeeze(cueBeam_xz(tx',wavenumber,x0,y0,z0,nx,ny,nz,dx,dy,dz));
tout=toc; raycount=nx*ny*nz*size(tx,1); rayspeed=raycount/tout; fprintf('%0.1f Mrays/s\n',rayspeed/1e6);
% convert image to decibel scale
img_xz_db=20*log10(img_xz./max(img_xz(:)));
figure(2);
imagesc([z0 z1],[x0 x1],img_xz_db,[-20 0]);
colorbar;
axis image
title('XZ cross-section: dB range');
xlabel('z-coordinate[m]'); ylabel('x-coordinate[m]');


