% finally! an equi-area projection from rectangular mesh to sphere and
% back.
% to be used as an accurate method of calculating total contrast 
mstruct=defaultm('eqaazim');
mstruct.angleunits='radians';
mstruct.origin=[pi/2 0 0];
%mstruct.geoid=[0.8 0];
mstruct=defaultm(mstruct);

density=1/16;
xv=[-1:density:1]*sqrt(2); yv=xv;
[x y]=meshgrid(xv,yv); x=x(:); y=y(:);
[lat2 lon2]=minvtran(mstruct,x,y);

xvis=cos(lon2).*cos(lat2);
yvis=sin(lon2).*cos(lat2);
zvis=sin(lat2);
idx=find(lat2>0);

clf;
subplot(1,2,1);
plot3(xvis(idx),yvis(idx),zvis(idx),'.'); axis equal
xlabel('-x-'); ylabel('-y-'); zlabel('-z-');
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
view(-20,70);
subplot(1,2,2);
plot(x(idx),y(idx),'.'); axis equal