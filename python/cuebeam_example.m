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
ImgResMultiplier=2;
% image resolution
dx = 1.0e-3/ImgResMultiplier;
dy = 1.0e-3/ImgResMultiplier;
dz = 1.0e-3/ImgResMultiplier;
% number of pixels in the imag
nx = 1; % ! Note, in this version, nx must be 1. 
ny = ceil(256*ImgResMultiplier);
nz = ceil(512*ImgResMultiplier); % note: in this example, the array extends along Z, and the depth is Y


% origin of the image
x0=1e-3; % should be slightly off-centre to avoid division by zero
y0=50.0e-3;
z0=-(nz/2)*dz;

% create an array
element_count = 16;
element_spacing = 4e-3;

array_elements_z=([1:element_count]*element_spacing); %#ok<NBRAK>
array_elements_z=array_elements_z-mean(array_elements_z); % centre around z-axis

% phase for each element - simple case. Simply bend the phases a bit.
test_elements_p=linspace(0.0,4,length(array_elements_z));
% phase for each element - use focussing. Overwrites previous phase
focus_x=0;
focus_y=0.2;
focus_z=-0.05;
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
field_py=cueBeam.beamsim_remote(pyargs('k',wavenumber,'elements_vectorized',elements_vectorized,'dy',dy,'dz',dz,'ny',uint32(ny),'nz',uint32(nz),'nx',uint32(nx),'z0',z0,'y0',y0));
field=ndarray2mat2(field_py);
t2=now;

% note that this is intended for the networked version
calculate_benchmark_figures=true;
if calculate_benchmark_figures
    % calculate performance stats
    ray_count=numel(field)*length(array_elements_z) ;
    t_roundtrip = (3600*24*(t2-t1));
    
    performance =  ray_count / t_roundtrip;
    datarate = 8*numel(field)/ t_roundtrip;
    fprintf('got %0.1f MRays/s, %0.3f MB/sec\n',performance*1e-6,datarate/1024/1024);
    fprintf('image size: %0.1f MB\n',length(uint8(field_py.base))  / 1024/1024);
    fprintf('turnover time: %0.1f sec\n',t_roundtrip);
end
%% display the field, dB scale
figure(1); clf;
y=linspace(y0,y0+dy*ny,size(field,1));
z=linspace(z0,z0+dz*nz,size(field,2));

field_decibels=20*log10(abs(field));
field_decibels=field_decibels-max(field_decibels(:));
handle_img=imagesc(z,y,field_decibels,[-40 0]); axis image
colorbar;
hold on;
plot(array_elements_z,zeros(size(array_elements_z)),'r.')
xlabel('z-axis,meters'); ylabel('y-axis,meters');
title('field amplitude, dB scale');
% draw a line at the focal point
hl=line([z0 z0+nz*dz],[focus_y focus_y]);
hl.Color='w'; hl.LineStyle = '-.';
%% display the field, line
figure(2); clf;
y_crossection=focus_y;
[~, y_crossection_idx]=min(abs(y-y_crossection));
y_dataline=field(y_crossection_idx,:);
y_dataline_abs=abs(y_dataline)-max(abs(y_dataline));
plot(z,y_dataline_abs)
xlabel('z-axis, meters'); ylabel('field amplitude at focal section, dB');
grid on
hold on;
% calculate the beam width and SNR
% find how many pixels are above -3dB?
PixelsAbove3dB=y_dataline_abs>-3;
% figure out the area 'nearby' by replicating the beam width to the left and right
PixelsNearbyR=circshift(PixelsAbove3dB,[0 1*sum(PixelsAbove3dB)]);
PixelsNearbyL=circshift(PixelsAbove3dB,[0 -1*sum(PixelsAbove3dB)]);
PixelsAwayR=circshift(PixelsAbove3dB,[0 2*sum(PixelsAbove3dB)]);
PixelsAwayL=circshift(PixelsAbove3dB,[0 -2*sum(PixelsAbove3dB)]);

PixelsAroundMain = PixelsAbove3dB+PixelsNearbyR+PixelsNearbyL+PixelsAwayR+PixelsAwayL;
plot(z,PixelsAroundMain)
BeamWidth=sum(PixelsAbove3dB)*dz;
EnergyRatio=sum(10.^(y_dataline_abs(PixelsNearby)/20))/sum(10.^(y_dataline_abs(PixelsAbove3dB)/20));
EnergyRatiodB=20*log10(EnergyRatio);

% make sure it is a zero-one vector only
PixelsAroundMain(PixelsAroundMain>0)=1; PixelsAroundMain(PixelsAroundMain<1)=0;
% mask off the main lobe signal and it's surrounding
y_dataline_abs_masked=y_dataline_abs;
y_dataline_abs_masked(PixelsAroundMain>0)=NaN;

[PeakSidelobeBeyondMain, PeakSidelobeBeyondMainLocIdx ] = max(y_dataline_abs_masked);
PeakSidelobeBeyondMainLocation=z(PeakSidelobeBeyondMainLocIdx);

% combine the effect of near-lobe energy and the peak off-axis lobe
SNR_Score = max(EnergyRatiodB,PeakSidelobeBeyondMain);

title(sprintf('beamwidth: %0.1f mm, SNR = %0.1f dB, worst sidelobe = %0.1f dB\n',BeamWidth*1e3,EnergyRatiodB,PeakSidelobeBeyondMain))


% check for grating lobes here 
% TODO: Check for grating lobes here

% plot over the areas considered
zLeft=z(find(PixelsAbove3dB,1,'first'));
zRight=z(find(PixelsAbove3dB,1,'last'));
hFillCore=fill([zLeft zLeft zRight zRight],[0 -30 -30 0],'g');
hFillCore.FaceAlpha = 0.1;
zLeftNearby=z(find(PixelsNearbyR,1,'first'));
zRightNearby=z(find(PixelsNearbyR,1,'last'));
hFillNearby=fill([zLeftNearby zLeftNearby zRightNearby zRightNearby],[0 -30 -30 0],'r');
hFillNearby.FaceAlpha = 0.1;
