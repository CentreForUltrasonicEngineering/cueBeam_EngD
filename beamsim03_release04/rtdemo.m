
function rtdemo

% real time demo for cudaBEAM code
% please note! 
% the images are normalized to a value found under the cursor

home;
fprintf('click on figures to move the beam around. Press r to perturb array element location . . . . \n')
fprintf('press z to increase sparsity and x to decrease sparsity\n')
in=load('hp_04_fatwing'); % load the array description
%in=load('hp_rect128'); % load the array description
hp_pos=in.hp_pos;
if isfield(in,'description')
    fprintf('%s\n',in.description);
end
% !!-------------  CHANGE SPARSITY HERE
sparse=1; % when sparse=1, array is dense at 3MHz
% !!------------- Change frequency range here
freq_tab=2.4e6:20e3:3.6e6;

% !!------------- Change Lambert map radius and resolution here
% resolution==density is in absolute distance between pixels, eg.0.5mm
r=15e-3; density=1e-3; 


perturbation_amount=0.1e-3; % that's for perturbing elements with 'r' key


fprintf('using %0.1f%% bandwidth and %d freq steps\n',100*(max(freq_tab)-min(freq_tab))/mean(freq_tab),length(freq_tab));
fprintf('\n               ');

section_zoom=4;  % obsolete
lambert_zoom=1;  % obsolete

hp_pos=hp_pos./2.*0.95e-3*sparse; % apply the "ensparsing"

fps_accumulator=0; fps_forget=0.2;


f1=figure(1); set(f1,'Interruptible','off')
set(f1,'WindowButtonDownFcn',@wbdcb)
set(f1,'KeyPressFcn',@keyPressHandler)
colormap(jet(256))
ah=subplot(2,1,1); set(ah,'drawmode','fast');
%ah = axes('DrawMode','fast');
vl=5600; frequency=3e6;

wavelength=vl/frequency;
k=single(2*pi/wavelength);
xbase0=-sqrt(2)+1e-7;

npts=ceil(6.283185307179586*r/density);
d=2*sqrt(2)/npts;    
n=ceil(2*sqrt(2)/d);
img_lambert=zeros(n,n,'single');

ahi=imagesc(img_lambert,[-20 0]); axis image;% colorbar;
title(sprintf('lambert equiareal azimuthal map at r=%0.1f mm',r*1e3));
% plot a line showing xz plane
hlinexz=line([0 n],[n/2+1 n/2+1],'Parent',ah,'color','w');

ahxz=subplot(2,1,2); set(ahxz,'drawmode','fast');


x0=single(-30e-3); y0=single(0); z0=single(0e-3); xe=single(60e-3); ye=single(0); ze=single(80e-3);

z0_peak_start=single(5e-3);
nx=uint32(ceil(xe/density)); ny=uint32(1); nz=uint32(ceil(ze/density));
nz_peak_start=ceil(z0_peak_start/density);
dx=single(density); dy=dx; dz=dx;

do_align=0; % saves on unused threads, but not neccesarly makes any useful work.
if do_align==1
align_x=24; align_z=16; % size of computing block
nx=uint32(align_x*ceil(nx/align_x));
nz=uint32(align_z*ceil(nz/align_z));
end
img_xz=zeros(nx,nz,'single'); 

ahi2=imagesc([z0 z0+ze],[x0 x0+xe],img_xz,[-20 0]); axis image;
title('xz plane - section'); xlabel('x[m]'); ylabel('z[m]')
% plot a line showing a shilouette of lambert sphere
a=(-pi/2:pi/16:pi/2); xss=r*cos(a); yss=r*sin(a);
hline_lambert=line(xss,yss,'Parent',ahxz,'color','w');
TPoints=hp_pos;
TPoints=[TPoints zeros(size(hp_pos,1),1)];
iiUpdateDisplay(0,0,r,1);

f2=figure(2); 
set(f2,'KeyPressFcn',@keyPressHandler)
set(f2,'Interruptible','off')
plot(hp_pos(:,1),hp_pos(:,2),'ro'); 
axis image;
xlim([-10e-3 10e-3]); ylim([-10e-3 10e-3]);
figure(1);   
fprintf('running benchmark first . . .  \n                \n\n');

for zz=(-pi/2):(pi/500):(pi/2)
   iiUpdateDisplay(sin(zz)*r,0,cos(zz)*r,0); drawnow; 
end


function keyPressHandler(src,event)
    if strcmp(event.Key,'r')
        % randomize array elem positions     
        hp_pos=hp_pos+perturbation_amount*randn(size(hp_pos));            
    end 
    if strcmp(event.Key,'z')
        hp_pos=hp_pos.*1.01;
    end
    if strcmp(event.Key,'x')
        hp_pos=hp_pos.*(1/1.01);
    end 
    
    
        TPoints=hp_pos;
        TPoints=[TPoints zeros(size(hp_pos,1),1)];
        figure(2); 
        plot(hp_pos(:,1),hp_pos(:,2),'ro'); 
         axis image;
        xlim([-10e-3 10e-3]); ylim([-10e-3 10e-3]);
        figure(1);
        
    wbmcb(src,event); % updates display
end

function iiUpdateDisplay(fx,fy,fz,peak_flag)        
    tic
           % focus the tx at lambert_xyz
           Distance = sqrt((TPoints(:,1)-fx).^2+(TPoints(:,2)-fy).^2+(TPoints(:,3)-fz).^2);  %!!NOTE: This steers per-transmiter point, not per element if there are more than 1 points per element. May need to adjust this later on!           
           SVect=exp(1i*2*pi*Distance/wavelength); % steering vector. need to update that in the callback too!
           tx=[hp_pos(:,1) hp_pos(:,2) zeros(size(hp_pos,1),1) zeros(size(hp_pos,1),1) abs(SVect(:)) angle(SVect(:))];
           img_xz=zeros(nx,nz,'single'); 
           img_lambert=zeros(n+1,n+1,'single');
           for f=freq_tab
               wavelength=vl/f; k=single(2*pi/wavelength);
               SVect=exp(1i*2*pi*Distance/wavelength); % steering vector. need to update that in the callback too!
               tx=[hp_pos(:,1) hp_pos(:,2) zeros(size(hp_pos,1),1) zeros(size(hp_pos,1),1) abs(SVect(:)) angle(SVect(:))];
           img_lambert=img_lambert+beamsim_lambert_cu(single(tx)',single(k),single(r),single(density));           
           img_xz=img_xz+squeeze(beamsim0302_xz_cu(single(tx)',k,x0,y0,z0,nx,ny,nz,dx,dy,dz));
           end
           if peak_flag==0; % peak-flag=0 - take from labmert, xz otherwise
                peak=max(img_lambert(:)); 
           else
              % peak=max(max(img_xz(:,nz_peak_start:end)));
              % new: normalize to the point under the cursor
              % find x,y of the fy,fz under cursor
              pt_y=ceil((fx-x0)./dx);
              pt_x=ceil((fz-z0)./dx);
              peak=img_xz(pt_y,pt_x);
           end
           
           img_labmert=20*log10(img_lambert./peak);
           img_xz=20*log10(img_xz./peak);
           %set(ah,'CLim',[-20 0]);
           
           set(ahi,'CData',img_labmert);                                           
           set(ahi2,'CData',img_xz);
           
           %imagesc(imgl,[-20 0]);
           %drawnow;
           %set(ah,'drawmode','fast');
           %xdat = [xinit,cp(1,1)];
           %ydat = [yinit,cp(1,2)];
           %set(hl,'XData',xdat,'YData',ydat);drawnow
           tout=toc;
           fps=1/tout;
           fps_accumulator_new=(1-fps_forget)*fps_accumulator+fps_forget*fps;
           if (abs(fps_accumulator_new-fps_accumulator)/fps_accumulator)>0.005 % 0.5% accuracy required
               fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\bmeasuringFPS ')
               fps_accumulator=fps_accumulator_new;
               fps_forget=0.99*fps_forget;
              % title(sprintf('fps forget rate %0.5f',fps_forget));
           else
               fps_accumulator=fps_accumulator_new;
               fps_forget=0.999*fps_forget;
               fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b   %05.1f FPS ',fps_accumulator);               
           end
end
  function wbdcb(src,evnt)
     if strcmp(get(src,'SelectionType'),'normal')
        set(src,'pointer','circle')
        cp = get(ah,'CurrentPoint');
        %xinit = cp(1,1);yinit = cp(1,2);
        %hl = line('XData',xinit,'YData',yinit,...
        %'Marker','p','color','b');
        set(src,'WindowButtonMotionFcn',@wbmcb)
        set(src,'WindowButtonUpFcn',@wbucb)
     end
  end   
        function wbmcb(src,evnt)
           tic
           
           cp = get(ah,'CurrentPoint');
           %fprintf('\b\b\b\b\b\b\b\b\b%04.0f:%04.0f',cp(1,1),cp(1,2))
           % calculate xbase, ybase
           ix=cp(1,2); iy=cp(1,1);
           if (ix>0&&ix<n&&iy>0&&iy<n) % means that lambert window is updating
            xbase=ix*d+xbase0;
            ybase=iy*d+xbase0;
            % transform into lambert sphere
            rho2=xbase*xbase+ybase*ybase;
            if (rho2>2)
                return
            end 
            rhoi=1/sqrt(rho2);
            cosl=-ybase*rhoi;
            cosphi=sqrt(rho2-rho2*rho2/4); 
            lambert_x=r*cosl*cosphi;
            sinl=xbase*rhoi;
            lambert_y=r*sinl*cosphi;
            sinphi=1-rho2/2;
            lambert_z=r*sinphi;           
            iiUpdateDisplay(lambert_x,lambert_y,lambert_z,0);
           end
          cp=get(ahxz,'CurrentPoint');
          ix=cp(1,2); iz=cp(1,1);
%           if (ix>0&&ix<nx&&iz>0&&iz<nz) % means that xz plane is updating
%               focus_y=0;
%               focus_x=ix*density+x0;
%               focus_z=iz*density+z0;
%               iiUpdateDisplay(focus_x,focus_y,focus_z,1);
%           end
          if (ix>x0&&ix<(x0+xe)&&iz>z0&&iz<(z0+ze)) % means that xz plane is updating
              focus_y=0;
              focus_x=ix;
              focus_z=iz;
              iiUpdateDisplay(focus_x,focus_y,focus_z,1);
          end

        end
   
        function wbucb(src,evnt)
           if strcmp(get(src,'SelectionType'),'alt')
              set(src,'Pointer','arrow')
              set(src,'WindowButtonMotionFcn','')
              set(src,'WindowButtonUpFcn','')
           else
              return
           end
        end
  end

