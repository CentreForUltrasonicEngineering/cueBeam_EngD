% cueBeam linear array
% This script demonstrates how to use the cueBEAM.
% Copyright Jerzy Dziewierz, University of Strathclyde, 2008-2015

% set up enviroment parameters. All units in SI
%probe.frequency=5e6;
%enviroment.wave_velocity=5600; % steel
figureOffset=74; % this is to avoid collisions with other software

% prepare output dir
try
    [status,message,messageid]=mkdir(simulation.prefix);
catch E
end

% enable line smoothing to make plots nicer
try
if simulation.LineSmoothing
    set(0,'DefaultLineLineSmoothing','on')
    set(0,'DefaultPatchLineSmoothing','on')
else
    set(0,'DefaultLineLineSmoothing','off')
    set(0,'DefaultPatchLineSmoothing','off')
end
catch E % silent fallover
end

enviroment.wavelength=enviroment.wave_velocity/probe.frequency;
enviroment.wavenumber=single(2*pi/enviroment.wavelength);
enviroment.dx_simulation=enviroment.wavelength/7; % points per wavelength for discrete rayleigh integral

% set XYZ of focal point
% beam.alpha_rotation=pi/12;
% beam.alpha_rotation=-0.1*0.5*pi;
% beam.alpha_rotation=0;
% beam.focal_distance=50e-3;

beam.focal_x=beam.focal_distance*sin(beam.alpha_rotation);
beam.focal_y=0;
beam.focal_z=beam.focal_distance*cos(beam.alpha_rotation);

%beam.display_limit_db=-20;


% calculate effective probe parameters
probe.g=probe.d-probe.e; % distance between elements
probe.A=probe.n*probe.d;

probe.active_nearfield=probe.A.^2*probe.frequency/4/enviroment.wave_velocity;
probe.passive_nearfield=probe.W.^2*probe.frequency/4/enviroment.wave_velocity;
probe.active_elementsize_lambda=probe.e/enviroment.wavelength;
probe.passive_elementsize_lambda=probe.W/enviroment.wavelength;
probe.active_aperture_lambda=probe.A/enviroment.wavelength;

% calculate active_aperture_focussing
bd=@(FL,v,f,A)(1.02*FL*v/f/A);
probe.active_aperture_focussing=bd(beam.focal_distance,enviroment.wave_velocity,probe.frequency,probe.A);

if simulation.verbose
    fprintf('-------: calculated characteristics :-------\n');
    fprintf('probe aperture: %0.1f x %0.1f mm\n',probe.A*1e3,probe.W*1e3);
    fprintf('probe active plane element size: %0.1f lambda\n',probe.active_elementsize_lambda);
    fprintf('probe passive plane element size: %0.1f lambda\n',probe.passive_elementsize_lambda);
    fprintf('probe active plane aperture: %0.1f lambda\n',probe.active_aperture_lambda);
    fprintf('probe active plane nearfield: %0.1f mm\n',probe.active_nearfield*1e3);
    fprintf('probe passive plane nearfield: %0.1f mm\n',probe.passive_nearfield*1e3);
    fprintf('probe active plane -6dB focussing power estimate: %0.1fmm\n',probe.active_aperture_focussing*1e3);
end % simulation.verbose

% note: passive plane has 'natural' focus at approximately the nearfield-farfield
% boundary.


% generate element locations
% note: this is a generator for element locations of a linear probe. To be
% precise, the probe centre points are created here.
probe.x=linspace(-probe.A/2,probe.A/2,probe.n);
probe.y=zeros(size(probe.x));
probe.z=zeros(size(probe.x));

% calculate delay laws to focus the probe at focal point
for idx=1:probe.n
    probe.distanceToToFocalPoint(idx)=sqrt((beam.focal_x-probe.x(idx)).^2+(beam.focal_y-probe.y(idx)).^2+(beam.focal_z-probe.z(idx)).^2);
    probe.ToF(idx)=probe.distanceToToFocalPoint(idx)/enviroment.wave_velocity;

    % You can implement Your own apodisation type here if you want
    if probe.apodisationtype==cueBeam.ApodisationType.None
        probe.apodisation(idx)=1;
    elseif probe.apodisationtype==cueBeam.ApodisationType.RaisedCosine
        probe.apodisation(idx)=probe.apodisationParameter1+(1-probe.apodisationParameter1)*sin((idx-1)/(probe.n-1)*pi).^probe.apodisationParameter2;
    else
        error('unimplemented apodisation type');
    end
    
end

% correct delay laws to be negative only
probe_max_ToF=max([probe.ToF]);
probe.ToF=probe.ToF-probe_max_ToF;

if simulation.doplots
    figure(figureOffset+5); clf;
    set(gcf,'Name','Firing delays','NumberTitle','off');
    bar([probe.ToF]); title('firing delays applied to array elements');
    ylabel('delay[s]'); xlabel('element number[-]');
    if simulation.doprints
        %print('-dpng',simulation.printresolution,sprintf('%s\\%sdelays.png',simulation.prefix,simulation.prefix));        
        cueBeam.myaa(simulation,sprintf('%s\\%sdelays.png',simulation.prefix,simulation.prefix));
    end
    % export delays to PZFlex format
    % this is for putting the delays into a PZFlex simulation
    try
        fout=fopen(sprintf('%s\\%s_pzflex_delays.in',simulation.prefix,simulation.prefix),'w+');
        for idx=1:probe.n
            fprintf(fout,'symb tshift%d = %e\n',idx,  probe.ToF(idx));
        end
        fclose(fout);
    catch E
        warning('problem exporting delays to PZFlex');
    end
end

% to simulate element directivity, use multiple points for each probe
% element. This is akin to numerical integration over a plane.
tx=[]; % reset emitter points description
for idx=1:probe.n
    % width direction
    npts_x=ceil(probe.e/enviroment.dx_simulation);
    px=linspace(probe.x(idx)-probe.e/2,probe.x(idx)+probe.e/2,npts_x); % x-coordinate points
    npts_y=ceil(probe.W/enviroment.dx_simulation);
    py=linspace(probe.y(idx)-probe.W/2,probe.y(idx)+probe.W/2,npts_y);
    [pxx pyy]=meshgrid(px,py); pxx=pxx(:); pyy=pyy(:);
    pzz=zeros(size(pyy));
    pzeros=zeros(size(pyy));
    element_steering_phase=2*pi*probe.distanceToToFocalPoint(idx)/enviroment.wavelength;
    pff=element_steering_phase+zeros(size(pxx));
    pfa=probe.apodisation(idx)+zeros(size(pxx));
    tx=[tx; pxx pyy pzz pzeros pfa pff];
end
tx=single(tx);
% save
probe.tx=tx;

%% make lambert map image

% define image parameters
lambert_radius=single(beam.focal_distance);
simulation.lambert_radius=lambert_radius;
lambert_map_density=single(simulation.lambert_map_density);


% !! Note: This is where the ray-tracing takes place.
tic
[img_lambert lambert_x lambert_y lambert_z]=cueBeam.cueBeam_lambert(tx',enviroment.wavenumber,lambert_radius,lambert_map_density);
tbenchmark=toc;
% ray tracing complete.

% uncomment for benchmark:
% raycount=numel(img_lambert)*size(tx,1); rayspeed=raycount/tout; fprintf('%0.1f Mrays/s\n',rayspeed/1e6);
% convert image to decibel scale
img_db=20*log10(img_lambert./max(img_lambert(:)));
beam.img_lambert=img_db;
beam.img_lambert_x=lambert_x;
beam.img_lambert_y=lambert_y;
beam.img_lambert_z=lambert_z;
lambert_coords=linspace(-90,90,size(img_db,1));
beam.img_lambert_coords=lambert_coords;

% display lambert map image

if simulation.doplots
    figure(figureOffset+1); clf;
    set(gcf,'Name','Lambert map, dB range','NumberTitle','off');
    clf;
    imagesc(lambert_coords,lambert_coords,img_db,[beam.display_limit_db 0]); colorbar; axis image;
    title(sprintf('Lambert map : dB range\npixel edge length: %0.2f mm',1e3*lambert_map_density));
    xlabel('angle[deg]'); ylabel('angle[deg]');
    if simulation.doprints
        %print('-dpng',simulation.printresolution,sprintf('%s\\%slambert_map.png',simulation.prefix,simulation.prefix));
        cueBeam.myaa(simulation,sprintf('%s\\%slambert_map.png',simulation.prefix,simulation.prefix));
    end
    
end
%%
% make XZ cross-section image
% define XZ image parameters
resolution=single(simulation.xy_resolution);

x0=single(-simulation.xy_x_extent); x1=single(simulation.xy_x_extent); dx=resolution; nx=uint32(ceil((x1-x0)/dx));
y0=single(0); y1=single(0); dy=resolution; ny=uint32(ceil((y1-y0)/dy)+1);
z0=single(simulation.xy_z0); z1=single(simulation.xy_z_extent); dz=resolution; nz=uint32(ceil((z1-z0)/dz));

% calculated image coordinates
xpoints=x0:dx:(x0+dx*single(nx-1));
ypoints=y0:dy:(y0+dy*single(ny-1));
zpoints=z0:dz:(z0+dz*single(nz-1));

% perform beam simulation
% !! Note: This is where the ray-tracing takes place.
tic;
img_xz=squeeze(cueBeam.cueBeam_xz(tx',enviroment.wavenumber,x0,y0,z0,nx,ny,nz,dx,dy,dz));
tbenchmark=toc;
% ray tracing complete.

% uncomment for benchmark
%raycount=nx*ny*nz*size(tx,1); rayspeed=raycount/tbenchmark; fprintf('%0.1f Mrays/s\n',rayspeed/1e6);

% convert image to decibel scale
img_xz_db=20*log10(img_xz./max(img_xz(:)));

beam.img_xz_xpoints=xpoints;
beam.img_xz_ypoints=ypoints;
beam.img_xz_zpoints=zpoints;
beam.beam_img_xz=img_xz_db;

if simulation.doplots
    figure(figureOffset+2); clf;
    set(gcf,'Name','XZ cross-section','NumberTitle','off');
    imagesc([z0 z1],[x0 x1],img_xz_db,[beam.display_limit_db 0]);
    colorbar;
    title('XZ cross-section: dB range');
    xlabel('z-coordinate[m]'); ylabel('x-coordinate[m]');
    hold on;
    plot(tx(:,3),tx(:,1),'go');
    axis image
    set(gca,'YDir','normal')
    if simulation.doprints
        %print('-dpng',simulation.printresolution,sprintf('%s\\%sbeam_xy.png',simulation.prefix,simulation.prefix));
        cueBeam.myaa(simulation,sprintf('%s\\%sbeam_xy.png',simulation.prefix,simulation.prefix));
    end
    % plot the same as contour
    figure(figureOffset+9); clf;
    set(gcf,'Name','XZ cross-section contour','NumberTitle','off');
    contours=[-20 -6 -3];
    [C h]=contour(zpoints,xpoints,img_xz_db,contours);
    set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
    set(h,'LineWidth',2);
    title('XZ cross-section, contour: dB range');
    xlabel('z-coordinate[m]'); ylabel('x-coordinate[m]');
    hold on;
    plot(tx(:,3),tx(:,1),'go');
    axis image;
    ylim([min(xpoints) max(xpoints)])
    grid on;
    colorbar;
    set(gca,'YDir','normal')
    if simulation.doprints
        %print('-dpng',simulation.printresolution,sprintf('%s\\%sbeam_xy_contour.png',simulation.prefix,simulation.prefix));
        cueBeam.myaa(simulation,sprintf('%s\\%sbeam_xy_contour.png',simulation.prefix,simulation.prefix));
    end
end
%% prepare 3D display of XZ

if simulation.do3Dplot1
    figure(figureOffset+3); clf;
    set(gcf,'Name','3D preview','NumberTitle','off');
    
    plot3(tx(:,1),tx(:,2),tx(:,3),'r.'); hold on;
    
    % create xz mesh
    [xpp zpp]=meshgrid(xpoints,zpoints); xpp=xpp(:); zpp=zpp(:); ypp=zeros(size(xpp));
    tri=delaunay(double(xpp),double(zpp));
    cpp=img_xz_db'; cpp=cpp(:);
    cpp_alpha=cpp; cpp_alpha(cpp>-20)=1; cpp_alpha(cpp<=-20)=0.1;
    cpp_alpha(1)=0;
    hs=trisurf(tri,xpp,ypp,zpp,cpp);
    set(hs,'LineStyle','none');
    set(hs,'FaceVertexAlphaData',cpp_alpha);
    set(hs,'FaceAlpha','interp')
    shading interp
    img_db_alpha=zeros(size(img_db));
    img_db_alpha(img_db>beam.display_limit_db)=1;
    img_db_alpha(img_db<beam.display_limit_db)=0.03;
    img_db_alpha(1,1)=0;
    hl=surf(lambert_x,lambert_y,lambert_z,double(img_db));
    set(hl,'LineStyle','none');
    set(hl,'AlphaData',img_db_alpha)
    set(hl,'FaceAlpha','interp')
    set(hl,'AlphaDataMapping','scaled')
    
    caxis([beam.display_limit_db 0])
    axis image;
    camproj persp
    grid on;
    if simulation.doprints
        %print('-dpng',simulation.printresolution,sprintf('%s\\%sperspective.png',simulation.prefix,simulation.prefix));
        cueBeam.myaa(simulation,sprintf('%s\\%sperspective.png',simulation.prefix,simulation.prefix));
    end
end % simulation.do3Dplot

%% do 3D beam and export it to Voreen

if simulation.do3DBeam
    if simulation.verbose
        fprintf('3D volume slices left: -----');
    end
    % extend y to same extents as x, leave z as was
    ypoints=xpoints;
    datacube=zeros(nx,nx,nz,'single');
    for idx_y=1:length(ypoints)
        if simulation.verbose
            fprintf('\b\b\b\b\b%05d',length(ypoints)-idx_y);
        end
        yo_local=single(ypoints(idx_y));
        img_xz=squeeze(cueBeam.cueBeam_xz(tx',enviroment.wavenumber,x0,yo_local,z0,nx,ny,nz,dx,dy,dz));
        datacube(idx_y,:,:)=img_xz;
    end
    % convert to log scale
    datacube=20*log10(abs(datacube).*(1/max(datacube(:))));
    
    if simulation.verbose
        fprintf('\b\b\b\b\bsaving..');
    end
    cueBeam.ExportToVoreen(simulation.prefix,sprintf('%sVolume',simulation.prefix),datacube,[beam.display_limit_db 0],[dx dz dz]);
    if simulation.verbose
        fprintf('\b\b\b\b\b\b\b\b');
        fprintf('\b\b\b\b\b\b');
        fprintf('completed, now load/refresh the file in voreen\n')
    end
end

%% beam section on XZ plot
if simulation.doXZBeamSection
    [location_error idx_z]=min(abs(simulation.XZBeamSection_Z-zpoints)); % find closest line that HAS been calculated
    XLine=beam.beam_img_xz(:,idx_z);
    % normalize line against itself
    XLine=XLine-max(XLine);
    [Xintersect Yintersect]=cueBeam.curveintersect([min(xpoints) max(xpoints)],[-3 -3],xpoints,XLine);
    
    
    
    if length(Xintersect)==2
        beamwidth=Xintersect(2)-Xintersect(1);
        str=sprintf('-3dB beam width (linear): %0.1f mm',beamwidth*1e3);
        beam.beamwidth_XZ=beamwidth;
        
        % get mean amplitude of everything EXCEPT the main lobe
        [location_error idx_1]=min(abs(Xintersect(1)-xpoints));
        [location_error idx_2]=min(abs(Xintersect(2)-xpoints));
        XLine_cut=10.^(XLine/20);
        XLine_cut(idx_1:idx_2)=[];
        IntegratedSidelobeXLine=20*log10(mean(XLine_cut));
        beam.IntegratedSidelobe_XLine=IntegratedSidelobeXLine;
        
        beam.peakSidelobe_XLine=cueBeam.GetPeakSidelobeValue(XLine);
        
        
    else
        str=sprintf('no -3dB intersection or multiple intersections');
        beam.beamwidth_XZ=NaN;
        beam.IntegratedSidelobe_XLine=NaN;
    end
    if simulation.doplots
        figure(figureOffset+6); clf;
        set(gcf,'Name','XZ plane section','NumberTitle','off')
        hpx=plot([min(xpoints) max(xpoints)],[-3 -3],'g',xpoints,XLine,'k',Xintersect,Yintersect,'ro'); grid on;
        set(hpx,'linewidth',2);
        
        ylim([beam.display_limit_db 1]);
        xlabel('x-dimension[m]'); ylabel('amplitude vs. peak[dB]');
        title(sprintf('XZ plane Z-section at z=%0.1f mm\n%s',zpoints(idx_z)*1e3,str));
        
        if simulation.doprints
            %print('-dpng',simulation.printresolution,sprintf('%s\\%sXZ_plane_section.png',simulation.prefix,simulation.prefix));
            cueBeam.myaa(simulation,sprintf('%s\\%sXZ_plane_section.png',simulation.prefix,simulation.prefix));
        end
        
    end
    if simulation.verbose
        fprintf('XZ plane section beam width: %0.1f mm\n',beam.beamwidth_XZ*1e3);
        fprintf('XZ plane section integrated sidelobe: %0.1f dB\n',beam.IntegratedSidelobe_XLine);
        fprintf('XZ plane section peak sidelobe: %0.1f dB\n',beam.peakSidelobe_XLine);
    end
end
%%
if simulation.doLambertSection
    [location_error idx_steer]=min(abs(0-beam.img_lambert_coords));
    SteeringPlaneLine=beam.img_lambert(idx_steer,:);
    % normalize line against itself
    SteeringPlaneLine=SteeringPlaneLine-max(SteeringPlaneLine);
    [Xintersect Yintersect]=cueBeam.curveintersect([-90 90],[-3 -3],beam.img_lambert_coords,SteeringPlaneLine);
    if length(Xintersect)==2
        beamwidth=Xintersect(2)-Xintersect(1);
        str=sprintf('-3dB beam width (angular): %0.1f degrees',beamwidth);
        beam.beamwidth_lambert=beamwidth;
        
        % get mean amplitude of everything EXCEPT the main lobe
        [location_error idx_1]=min(abs(Xintersect(1)-beam.img_lambert_coords));
        [location_error idx_2]=min(abs(Xintersect(2)-beam.img_lambert_coords));
        SteeringPlaneLine_cut=10.^(SteeringPlaneLine/20);
        SteeringPlaneLine_cut(idx_1:idx_2)=[];
        IntegratedSidelobe=20*log10(mean(SteeringPlaneLine_cut));
        beam.IntegratedSidelobe_lambert=IntegratedSidelobe;
        
        beam.peakSidelobe_lambert=cueBeam.GetPeakSidelobeValue(SteeringPlaneLine);
        
    else
        str=sprintf('no -3dB intersection or multiple intersections');
        beam.beamwidth_lambert=NaN;
        beam.IntegratedSidelobe_lambert=NaN;
    end
    
    if simulation.doplots
        figure(figureOffset+7); clf;
        set(gcf,'Name','Lambert sphere section','NumberTitle','off')
        hpl=plot([-90 90],[-3 -3],'g',beam.img_lambert_coords,SteeringPlaneLine,'k',Xintersect,Yintersect,'ro'); grid on;
        set(hpl,'linewidth',2);
        ylim([beam.display_limit_db 1]);
        xlabel('steering angle[degrees]'); ylabel('amplitude vs. peak[dB]');
        
        title(sprintf('Lambert sphere section in steering plane\nsphere radius %0.1f mm\n%s',simulation.lambert_radius*1e3,str));
        
        if simulation.doprints
            %print('-dpng',simulation.printresolution,sprintf('%s\\%sXZ_Lambert_section.png',simulation.prefix,simulation.prefix));
            cueBeam.myaa(simulation,sprintf('%s\\%sXZ_Lambert_section.png',simulation.prefix,simulation.prefix));
        end
        
    end
    
    if simulation.verbose
        fprintf('Lambert sphere section beam width: %0.1f degrees\n',beam.beamwidth_lambert);
        fprintf('Lambert sphere section integrated sidelobe: %0.1f dB\n',beam.IntegratedSidelobe_lambert);
        fprintf('Lambert sphere section peak sidelobe: %0.1f dB\n',beam.peakSidelobe_lambert);
    end
end
