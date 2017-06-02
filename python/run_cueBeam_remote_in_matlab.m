% Demonstration: how to use cueBeam python module from Matlab


% note: if the python source code is modified, it has to be explicitly re-loaded in Matlab
% one has to do 'clear classes' in such case.
% clear classes
cueBeam=py.importlib.import_module('cueBeamCore3');
py.importlib.reload(cueBeam);
%% create an example scenario (world)
frequency = 1e6; %Hz
wave_velocity_in_medium = 1450; % m/s; water=1450, plastic = 2800, steel = 5600

wavelength=wave_velocity_in_medium/frequency; 
wavenumber=1/wavelength;
%%
ImgResMultiplier=0.25;
% image resolution
dx = 1.0e-3/ImgResMultiplier; 
dy = 1.0e-3/ImgResMultiplier;
dz = 1.0e-3/ImgResMultiplier;
% number of pixels in the image
nx = 1; 
ny = 512*ImgResMultiplier;
nz = 512*ImgResMultiplier; % note: in this example, the array extends along Z, and the depth is Y

% origin of the image
x0=0.1e-3; % should be slightly off-centre to avoid division by zero
y0=5.0e-3;
z0 = -(dz/2)*nz;

% create an array
element_count = 32;
element_spacing = 3e-3;

array_elements_z=([1:element_count]*element_spacing); %#ok<NBRAK>
array_elements_z=array_elements_z-mean(array_elements_z); % centre around z-axis

% phase for each element - simple case. Simply bend the phases a bit.
test_elements_p=linspace(0.0,4,length(array_elements_z));
% phase for each element - use focussing.
focus_x=0;
focus_y=0.3;
focus_z=0.2;
for idxE=1:length(array_elements_z)    
    distance_element_to_focus=sqrt(((focus_x-0).^2)+((focus_y-0).^2)+((focus_z-array_elements_z(idxE))^2));
    phase_distance=distance_element_to_focus/wavelength;
    test_elements_p(idxE)=phase_distance;    
end

% encode the elements into an element description vector
elements_vectorized=[]; %elements_vectorized=[x,y,z,amplitude,phase,zero];
for idxE=1:length(array_elements_z)    
    x=0; y=0; z=array_elements_z(idxE); amp=1; phase=test_elements_p(idxE); tmp=0;
    elements_vectorized=[elements_vectorized x y z amp phase tmp]; %#ok<AGROW>
end

% call beamsim to do the work 
t1=now;
field=ndarray2mat(cueBeam.beamsim_remote(pyargs('k',wavenumber,'elements_vectorized',elements_vectorized,'dy',dy,'dz',dz,'ny',uint32(ny),'nz',uint32(nz),'z0',z0)));
t2=now;

% calculate performance stats
ray_count=numel(field)*length(array_elements_z) ;
t_roundtrip = (3600*24*(t2-t1));

performance =  ray_count / t_roundtrip;
datarate = 8*numel(field)/ t_roundtrip;
fprintf('got %0.1f MRays/s, %0.3f MB/sec\n',performance*1e-6,datarate/1024/1024);

%% display the field, dB scale
figure(1); clf;
y=linspace(y0,y0+dy*ny,size(field,1));
z=linspace(z0,z0+dz*nz,size(field,2));
field_decibels=20*log10(abs(field));
field_decibels=field_decibels-max(field_decibels(:)); 
handle_img=imagesc(z,y,field_decibels,[-40 0]); axis image
hold on;
plot(array_elements_z,zeros(size(array_elements_z)),'r.')
xlabel('z-axis,meters'); ylabel('y-axis,meters');


