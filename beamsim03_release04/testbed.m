% test bench for the beamsim03Cuda
% 16x16 2d array lying on z=0 plane, focused at 15 x=y=0 z=15 wavelengths
% do everything with single precision

%% define test parameters: 
% density controls image resolution. The lower the density, the higher the
% resolution and path count

density=1/550;

% please note: my GT9500 has theoretical capability of 134.4 GFLOPs, 
% the GTX480, in addition to power of 1344.96GFlops has some very important
% hardware features such as on-chip memory cache and larger register file.
% This means that this code is expected to run more than 10x faster on a new
% graphics card.

%% set up parameters. Feel free to manipulate these.

freq=1; vl=1; % lambda=1;
wavelength=vl/freq;
xfocus=-1; yfocus=0; zfocus=6;
% define transducer
% this is a dense matrix square-packed 2D array, (0,0) centered
probe_size_x=16;
probe_size_y=16;
tx_x=[1:probe_size_x]*wavelength/2; tx_x=tx_x-mean(tx_x);
tx_y=[1:probe_size_y]*wavelength/2; tx_y=tx_y-mean(tx_y);
[tx_x tx_y]=meshgrid(tx_x,tx_y); tx_x=tx_x(:); tx_y=tx_y(:);
tx_z=zeros(size(tx_x));
tx_size=0.5*ones(size(tx_x)); % note that tx_size is not used for calculations for now due to dispute on exact algorithm
                              % in other words, the code assumes that
                              % element is ominidirectional.

% randomize element positions for fun!
tx_x=tx_x+randn(size(tx_x))*wavelength;
tx_y=tx_y+randn(size(tx_y))*wavelength;
%% steering vector
% modify only if you know what you are doing
Distance = sqrt((tx_x-xfocus).^2+(tx_y-yfocus).^2+(tx_z-zfocus).^2);  %!!NOTE: This steers per-transmiter point, not per element if there are more than 1 points per element. May need to adjust this later on!
SVect=exp(1i*2*pi*Distance/wavelength); % steering vector

%% assemble transmiter vector
tx=single([tx_x tx_y tx_z tx_size abs(SVect) angle(SVect)])';
%% image parameters : testbed2: z-plane at z=15 extending 7
img_density=single(density);
img_xsize=single(7.0);
img_ysize=single(7.0);
img_zsize=single(0);

x0=-img_xsize/2; 
y0=-img_ysize/2;
z0=single(zfocus); 

nx=uint32(ceil(img_xsize/img_density));
ny=uint32(ceil(img_ysize/img_density));
nz=uint32(1);

dx=img_density;
dy=img_density;
dz=img_density;

% align no. of pixels to GPU's preferable size. saves wasted processing
% power on the unused boundaries.
do_align=1; % update: it appears that this yields measurable, but neglible difference. it may be just as well disabled            
if do_align==1
align_x=16; align_y=24; % size of computing block
nx=uint32(align_x*ceil(nx/align_x));
ny=uint32(align_y*ceil(ny/align_y));
end

k=single(2*pi/wavelength);

% calculate no. of paths
paths_no=nx*ny*nz*length(tx_x);
home;

fprintf('image size: %d x %d px; %0.1f M paths\n',nx,ny,paths_no/1e6); 

%% compile sources - uncomment this as neccesary
% mex -g beamsim03_c.c % -g for debugging
% mex -O beamsim03_c.c % -O enables optimizations for production - approx. 15% faster
% nvc beamsim03_cup5
%% call the calculation procedure
fprintf('please be patient while i conduct benchmark . . .\n')
if 1==0 % do or don't do m-version
fprintf('--matlab--...'); 
tic; img_m=beamsim03_m(tx,k,x0,y0,z0,nx,ny,nz,dx,dy,dz); t_m=toc;
%t_m=50.89; % for density=1/240
fprintf('% 8.2f M paths/s . . . done\n',1e-6*double(paths_no)/t_m);
% --matlab--...    4.99 M paths/s . . . don
end
if 1==0 % do or don't c-mex version
fprintf('--c-------...');
tic; img_c=beamsim03_c(tx,k,x0,y0,z0,nx,ny,nz,dx,dy,dz); t_c=toc;
%t_c=26.78 % for density=1/240
fprintf('% 8.2f M paths/s . . . done\n',1e-6*double(paths_no)/t_c);
%--c-------...   15.95 M paths/s . . . done on a core i7 @ 2.8GHz
end

%dif=mean(img_m(:)-img_c(:)) 
%nt=t_m/t_c; fprintf('speedup=%0.1f; diff=%0.3e /pixel\n',nt,dif);
%img_c_l=20*log10(img_c./max(img_c(:)));

%%
fprintf('--cuda----...');
% warmup - only really needed on the first call of the new kernel. This
% compiles the new intermediate code to the binary code. Not needed in production code as your graphics card doesn't ever change.
img_cu=beamsim03_cup6(tx,k,x0,y0,z0,nx,ny,nz,dx,dy,dz);
n=64; % run cuda n times to measure the time better
tic;
for nc=1:n
    img_cu=beamsim03_cup6(tx,k,x0,y0,z0,nx,ny,nz,dx,dy,dz);
end
tt=toc;
t_cu=tt/n;
fprintf('% 8.2f M paths/s . . . done\n',1e-6*double(paths_no)/t_cu);
% --cuda----... 14943.00 M paths/s . . . done on a GTX480 without sincos()
% --cuda----... 14924.00 M paths/s . . . done on a GTX480 with sincos()
% --cuda----... 14945.00 M paths/s . . . with unroll 1
% --cuda----... 42611.00 M paths/s . . . on Nvidia GTX Titan, single GPU

%%
% verify result is the same as in matlab
fprintf('verification ');
difference=max(abs(img_m(:)-img_cu(:)))/max(img_cu(:));
if difference<1e-4 % the allowance for error is coarse in this kind of algorithm. 
                   % It's the metrics (spot size, sidelobe amplitude) that
                   % matter - i'll code them to be inside CUDA kernel too. 
                   % in particular, low accuracy device intristic sin and cos has been used,
                   % can be changed for higher accuracy later on at mnimal
                   % expense of speed
    fprintf('passed : L1-norm %0.1e (%0.1f dB) \n',difference,20*log10(difference));    
else
    warning('fail : L1-norm %0.1e',difference);
end 
% boast 
fprintf('speedup cuda/C %0.2fx\n',t_c/t_cu);
fprintf('speedup cuda/matlab %0.2fx\n',t_m/t_cu);
fprintf('equivalent RT FPS: %0.1f\n',1/t_cu);
%return
%% display result
if 1==1
figure(1); clf;
subplot(2,2,1)
img2=squeeze(img_m); img_log=20*log10(abs(img2)./max(abs(img2(:))));
xvect=(x0:dx:(x0+single(nx)*dx)); yvect=(y0:dy:(y0+single(ny)*dy)); zvect=(z0:dz:(z0+single(nz)*dz)); 
[xvm yvm]=meshgrid(xvect,yvect); xvm=xvm(:); yvm=yvm(:); zvm=zeros(size(xvm));
imagesc(yvect,xvect,img_log,[-20 0]); axis image; title(sprintf('Matlab: %0.2f[s]',t_m))

subplot(2,2,2);
img2=squeeze(img_c); img_log=20*log10(abs(img2)./max(abs(img2(:))));
xvect=(x0:dx:(x0+single(nx)*dx)); yvect=(y0:dy:(y0+single(ny)*dy)); zvect=(z0:dz:(z0+single(nz)*dz)); 
[xvm yvm]=meshgrid(xvect,yvect); xvm=xvm(:); yvm=zvm(:); zvm=zeros(size(xvm));
imagesc(yvect,xvect,img_log,[-20 0]); axis image;
title(sprintf('C MEX: %0.2f[s]',t_c));

subplot(2,2,3);
img2=squeeze(img_cu); img_log=20*log10(abs(img2)./max(abs(img2(:))));
xvect=(x0:dx:(x0+single(nx)*dx)); yvect=(y0:dy:(y0+single(ny)*dy)); zvect=(z0:dz:(z0+single(nz)*dz)); 
[xvm yvm]=meshgrid(xvect,yvect); xvm=xvm(:); yvm=zvm(:); zvm=zeros(size(xvm));
h=imagesc(yvect,xvect,img_log,[-20 0]); axis image;
title(sprintf('CUDA MEX: %0.2f[s]',t_cu));

end
fprintf('all done.\n');