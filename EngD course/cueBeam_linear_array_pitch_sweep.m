% cueBeam_linear_array_sweep
% demonstrates how to conduct a parameter sweep on one of the probe
% properties
% Copyright Jerzy Dziewierz, University of Strathclyde, 2008-2013

% note, all units SI system

% set up specimen
enviroment.wave_velocity=5600; % m/s, steel

% set up beam parameters
beam.alpha_rotation=20*pi/180; % radians
beam.focal_distance=50e-3; % mm
beam.display_limit_db=-20; % for display purposes only

% set up probe element locations.
probe=[]; % reset probe description
probe.frequency=5e6; % Hz, all other
probe.n=15; % number of elements
%probe.d=0.7e-3; % element pitch % note: this goes into parameter sweep
%probe.e=probe.d-0.1e-3; % element width
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

simulation.verbose=true;
simulation.doplots=false;
simulation.doprints=true; % requires  simulation.doplots=true, but disable to make script go faster
simulation.do3Dplot1=false; % note, might be slow on some computers
simulation.do3DBeam=false; % export 3D beam to Voreen

simulation.doXZBeamSection=false;
simulation.XZBeamSection_Z=beam.focal_distance;

simulation.doLambertSection=true;

simulation.prefix='probe1_';
simulation.printresolution='-r200';
% launch the simulation script
result=[];
sweep_pitch=0.2e-3:0.05e-3:1.5e-3;
for idx_sweep=1:length(sweep_pitch)
    simulation.verbose=false;
    
    probe.d=sweep_pitch(idx_sweep); % element pitch
    probe.e=probe.d-0.1e-3; % element width
    
    fprintf('step %d - %2.3f mm\n',idx_sweep,probe.d*1e3);
    simulation.prefix=sprintf('probe_pitchsweep_%2.3f_',probe.d*1e3);
    
    cueBeam.process_linear_array;
    
    result(idx_sweep).enviroment=enviroment;
    result(idx_sweep).probe=probe;
    result(idx_sweep).beam=beam;
    result(idx_sweep).simulation=simulation;
    
    
end

% save result

save('sweep_result','result');

%% display beam width and average sidelobe vs element pitch
bwl=[]; % beamwidth buffer
bsl=[]; % side lobe integrated
bpl=[]; % side lobe peak
for idx_sweep=1:length(sweep_pitch)
    bwl(idx_sweep)=result(idx_sweep).beam.beamwidth_lambert;
    bsl(idx_sweep)=result(idx_sweep).beam.IntegratedSidelobe_lambert;
    bpl(idx_sweep)=result(idx_sweep).beam.peakSidelobe_lambert;
end
figure(8); clf;
ax1=axes;
p1=plot(sweep_pitch,bwl,'bx-'); xlabel('element pitch[m]'); ylabel('beamwidth[degrees]');
set(p1,'linewidth',4);
set(p1,'parent',ax1);
set(ax1,'YColor','b');
set(ax1,'YLim',[0 30])
set(ax1,'XLim',[0 max(sweep_pitch)])

ax2=axes; 
p2=plot(sweep_pitch,bsl,'ro-',sweep_pitch,bpl,'go-','parent',ax2); 
set(p2,'linewidth',4);
set(ax2,'YLim',[-30 -10]);
set(ax2,'XLim',[0 max(sweep_pitch)]);
set(ax2,'Color','none','YAxisLocation','right');
set(ax2,'YColor','r');
set(ax2,'Box','off')
ax2l=ylabel('integrated sidelobe level[dB]');
set(ax2l,'parent',ax2);
grid on

title(sprintf('beamwidth and side lobe vs element pitch for steering of %0.1f degrees',beam.alpha_rotation*180/pi));


