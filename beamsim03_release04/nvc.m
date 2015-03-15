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

options='-gencode=arch=compute_10,code=sm_10 -gencode=arch=compute_10,code=compute_10 -gencode=arch=compute_20,code=sm_20 -gencode=arch=compute_20,code=compute_20';
options=[options ' --use_fast_math'];
%options=[options ' -keep'];
txt=sprintf('nvcc %s.cu %s -c -lcufft -lcudart -lcuda --ptxas-options=-v -Ic:\\MATLAB2010a\\extern\\include\\',filename,options);
system(txt)

% now compile into mexw64 / mexw32 
%mex_options='-g'; % to include debug info 
mex_options='-O'; % enable optimisation

n=getenv('CUDA_LIB_PATH');
mex(['-L' n],mex_options,'-lcudart','-lcufft','-lcuda',sprintf('%s.obj',filename));

% clean up
delete(sprintf('%s.obj',filename));