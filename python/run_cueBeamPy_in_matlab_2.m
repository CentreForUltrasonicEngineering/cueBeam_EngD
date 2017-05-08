clear classes
%% testing the contained version
%T:\git\cueBeam\python
cueBeam=py.importlib.import_module('cueBeamCore');
py.importlib.reload(cueBeam);
cs=cueBeam.CueBeamSolver;
wavelength=1e-3;
k=1/wavelength;
dy=1.0e-3;
dz=1.0e-3;
test_elements_z=([1:8]*3e-3); %#ok<NBRAK>
test_elements_z=test_elements_z-mean(test_elements_z);
test_elements_p=linspace(0.0,0.0,length(test_elements_z));
elements_vectorized=[]; %elements_vectorized=[1,2,3,4,5,6,7,8,9,10,11,12];
for idxE=1:length(test_elements_z)    
    x=0; y=0; z=test_elements_z(idxE); amp=1; phase=test_elements_p(idxE); tmp=0;
    elements_vectorized=[elements_vectorized x y z amp phase tmp]; %#ok<AGROW>
end

%field=cell2mat(cell(cs.beamsim_simpler(cs,pyargs('k',k,'elements_vectorized',elements_vectorized))))
q=cs.beamsim_simpler(cs,pyargs('y0',2e-3,'k',k,'ny',uint32(120),'nz',uint32(80),'elements_vectorized',elements_vectorized));
%%
qf=20*log10(abs(ndarray2mat(q)));
qf=qf-max(qf(:)); 
himg=imagesc(qf,[-40 0]); axis image

