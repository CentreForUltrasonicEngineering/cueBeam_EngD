% nvc(filename) compile cuda via intermediate .obj file
% '.cu' extension is already added
% 
% in case if nvcc complains about some .h or .lim files, simply add them to
% the "nvcc.profile" file where they belong, do not use any shady nvmex.pl 
% file from the forums (may contain execucable code!).
% 
% regarding the options, simply comment or uncomment relevant line as
% needed

% 2010-07 Jerzy Dziewierz, CUE, Strathclyde University
% Public domain, but please keep this short comment above

function nvc(filename)
% make cuda .obj file first

clear mex
options='-gencode=arch=compute_10,code=sm_10 -gencode=arch=compute_10,code=compute_10 -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20 -gencode=arch=compute_35,code=sm_35';
options=[options ' --use_fast_math'];

debug=0;
if debug
    mex_options='-g'; % to include debug info 
    options=[options ' -G'];
    debugflag='-Xcompiler "/Zi"';
    %options=[options ' -keep']; % note, this sprinks a lot of darnduff
else
    mex_options='-O'; % enable optimisation
    options=[options ' -O'];
    debugflag='';
end



txt=sprintf('nvcc %s.cu %s -c -lcufft -lcudart -lcuda --ptxas-options=-v %s -Xcompiler "/openmp"',filename,options,debugflag);
system(txt)

% now compile into mexw64 / mexw32 
n=getenv('CUDA_LIB_PATH');
mex(['-L' n],mex_options,'-lcudart','-lcufft','-lcuda',sprintf('%s.obj',filename));

% clean up
delete(sprintf('%s.obj',filename));