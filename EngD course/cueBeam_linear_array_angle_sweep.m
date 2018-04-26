% cueBeam_linear_array_angle_sweep
% demonstrates how to conduct a parameter sweep on one of the probe
% properties
% Copyright Jerzy Dziewierz, University of Strathclyde, 2008-2013

% note, all units SI system

% set up specimen
enviroment.wave_velocity=5600; % m/s, steel

% set up beam parameters
beam.alpha_rotation=20*pi/180; % radians % Note: This goes into sweep
beam.focal_distance=50e-3; % mm
beam.display_limit_db=-20; % for display purposes only

% set up probe element locations.
probe=[]; % reset probe description
probe.frequency=5e6; % Hz, all other
probe.n=15; % number of elements
probe.d=0.7e-3; % element pitch 
probe.e=probe.d-0.1e-3; % element width
probe.W=15e-3;
probe.apodisationtype=cueBeam.ApodisationType.None;
probe.apodisationParameter1=0.5; % used with cueBeamApodisation.RaisedCosine
probe.apodisationParameter2=0.5;

% set simulation options
simulation.lambert_map_density=1e-3;
simulation.xy_resolution=0.5e-3;
simulation.xy_z_extent=100e-3;
simulation.xy_z0=15e-3;% note: do not start z0 from zero,
% because the decibel range plots will be overwhelmed by amplitude of the sound close to the elements.
simulation.xy_x_extent=50e-3;

simulation.verbose=false;
simulation.doplots=true;
simulation.doprints=true; % exports figures for report, requires  simulation.doplots=true, disable to make script go faster
simulation.papersize=15; %cm
simulation.aa_figures=false;

% note, might be slow on some computers
simulation.do3Dplot1=false; % create a 3D view
simulation.do3DBeam=false; % export 3D beam shape to Voreen

simulation.doXZBeamSection=false;
simulation.XZBeamSection_Z=beam.focal_distance;

simulation.doLambertSection=true;

simulation.prefix='probe1_';
% launch the simulation script
result=[];
azimuths=-45:5:45;
fname=@(x)sprintf('probe_azimuthal_sweep_%2.0f_',x);
for idx_sweep=1:length(azimuths)
    simulation.verbose=false;
       
    
    simulation.prefix=fname(azimuths(idx_sweep));
    fprintf('step %d - %s\n',idx_sweep,simulation.prefix);
    beam.alpha_rotation = azimuths(idx_sweep) * pi/180;
    % launch the simulation script
    cueBeam.process_linear_array;
    
    % save updated data structures
    result(idx_sweep).enviroment=enviroment;
    result(idx_sweep).probe=probe;
    result(idx_sweep).beam=beam;
    result(idx_sweep).simulation=simulation;
end

% save result

save('sweep_result','result');

%% display beam width and average sidelobe vs azimuth
bwl=[];
bsl=[];
for idx_sweep=1:length(azimuths)
    bwl(idx_sweep)=result(idx_sweep).beam.beamwidth_lambert;
    bsl(idx_sweep)=result(idx_sweep).beam.IntegratedSidelobe_lambert;
end
figure(1); clf;
ax1=axes;
p1=plot(azimuths,bwl,'bx-'); xlabel('eazimuth[degrees]'); ylabel('beamwidth[degrees]');
set(p1,'linewidth',4);
set(p1,'parent',ax1);
set(ax1,'YColor','b');
set(ax1,'YLim',[0 30])
set(ax1,'XLim',[min(azimuths) max(azimuths)])

ax2=axes;
p2=plot(azimuths,bsl,'ro-','parent',ax2);
set(p2,'linewidth',4);
set(ax2,'YLim',[-30 -10]);
set(ax2,'XLim',[min(azimuths) max(azimuths)]);
set(ax2,'Color','none','YAxisLocation','right');
set(ax2,'YColor','r');
set(ax2,'Box','off')
ax2l=ylabel('integrated sidelobe level[dB]');
set(ax2l,'parent',ax2);
grid on

title(sprintf('beamwidth and side lobe vs azimuthal scan'));

%% png2avi
if simulation.doprints==true
    fname2=@(x)sprintf('%s\\%sbeam_xy.png',fname(x),fname(x));
    avi = VideoWriter( 'animation.avi');
    avi.FrameRate=5;
    open(avi);
    frameorder=[1:length(azimuths)];
    frameorder=[frameorder fliplr(frameorder)];
    frameorder=[frameorder fliplr(frameorder)];
    for idx_sweep=frameorder
        f1=fname2(azimuths(idx_sweep));
        fprintf('rendering %s...\n',f1);
        frame=imread(f1);
        writeVideo(avi,frame)
    end
    close(avi);
    % have system run default action on the file
    system('animation.avi');
end


