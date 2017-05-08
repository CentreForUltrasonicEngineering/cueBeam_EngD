% Do computation by calling python
%clear classes
cueBeam=py.importlib.import_module('cueBeamCore');
py.importlib.reload(cueBeam);
cs=cueBeam.CueBeamSolver;
cs.wavenumber=1/2e-3;
cs.elements.clear
element_z=(0:5)*5e-3;
element_z=element_z-mean(element_z);
nel=cs.TxElement;
for idx=1:length(element_z)            
    cs.elements.append(cueBeam.copy.copy(nel));    
    %cs.elements{idx}.setXYZ(pyargs('self',cs,'x',0,'y',0,'z',element_z(idx)));
    py.setattr(cs.elements{idx},'z',element_z(idx))
end

%%
tic;
cs.beamsim(cs);
t1=toc;
field = ndarray2mat(cs.rxPlane.pressurefield);
tout=toc;
raycount=double(cs.get_ray_count(cs));
raysPerSecond=raycount/tout;
fprintf('%0.2f MRays/s\n',raysPerSecond/1e6);
figure(1);  clf;
imagesc(abs(field)); axis image; set(gca,'YDir','normal');

%% testing the contained version
%T:\git\cueBeam\python
cueBeam=py.importlib.import_module('cueBeamCore');
py.importlib.reload(cueBeam);
cs=cueBeam.CueBeamSolver;
k=1/1e-3;
dy=1.0e-3;
dz=1.0e-3;
elements_vectorized=[1,2,3,4,5,6,7,8,9,10,11,12];

%field=cell2mat(cell(cs.beamsim_simpler(cs,pyargs('k',k,'elements_vectorized',elements_vectorized))))
q=cs.beamsim_simpler(cs,pyargs('elements_vectorized',elements_vectorized))

