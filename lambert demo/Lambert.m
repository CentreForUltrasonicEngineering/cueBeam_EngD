% lambert azimuthal equi-area sampling for fast beam characteristics
% determination
radius=50e-3;

req_resolution=12e-3; 
n_points_required=ceil(2*pi*radius/req_resolution)


density=2*sqrt(2)/n_points_required;

xv=[-1:density:(1)]*sqrt(2)+100*eps; yv=xv;
[x y]=meshgrid(xv,yv); x=x(:); y=y(:);
% if any(find(abs(xv)<0.0001))
%     fprintf('zero in xv\n')
% end
% lambert inverse projection: given x,y, on the map, what are the spherical
% coordinates?
% http://mathworld.wolfram.com/LambertAzimuthalEqual-AreaProjection.html

rho=sqrt(x.*x+y.*y); c=2*asin(0.5*rho);

% keep only c where imag(c)~0
idx=find(abs(imag(c))<eps);
length(idx)-length(x);
rho=rho(idx); c=c(idx); 
x=x(idx); y=y(idx);

phi1=pi/2; l0=0;
% phi is the standard parallel, l is the central longitude
phi=asin(cos(c).*sin(phi1)+(y.*sin(c).*cos(phi1))./rho);
l=l0+atan2((x.*sin(c)),(rho.*cos(phi1).*cos(c)-y.*sin(phi1).*sin(c)));

% visualize (phi,l) in cartesian system again
xvis=radius*cos(l).*cos(phi);
yvis=radius*sin(l).*cos(phi);
zvis=radius*sin(phi);
% use only positive zvis
idxp=find(zvis>0);
clf; subplot(1,2,1);
plot3(xvis(idxp),yvis(idxp),zvis(idxp),'.'); axis equal
xlabel('-x-'); ylabel('-y-'); zlabel('-z-');
xlim([-1 1]*2*radius); ylim([-1 1]*2*radius); zlim([-1 1]*2*radius);
view(-20,70);
% visualize which of the points are used in unwrapped mesh
subplot(1,2,2);
plot(x(idxp),y(idxp),'ro'); axis equal

return
plot(l); ylim([-pi pi]);
%%
mstruct=defaultm('eqaazim');
mstruct.angleunits='radians';
mstruct.origin=[pi/2 0 0];
%mstruct.geoid=[0.8 0];
mstruct=defaultm(mstruct);

density=1/9;
xv=[-1:density:1]*0.87056; yv=xv;
[x y]=meshgrid(xv,yv); x=x(:); y=y(:);
[lat2 lon2]=minvtran(mstruct,x,y);

xvis=cos(lon2).*cos(lat2);
yvis=sin(lon2).*cos(lat2);
zvis=sin(lat2);

clf;
plot3(xvis,yvis,zvis,'.'); axis equal
xlabel('-x-'); ylabel('-y-'); zlabel('-z-');
xlim([-1 1]); ylim([-1 1]); zlim([-1 1]);
view(-20,70);

% function [lat, lon] = eqaazimInv(mstruct, x, y)
% 
% [radius, D] = deriveParameters(mstruct);
% 
% rho = hypot(x/D, D*y);
% 
% indx = find(rho ~= 0);
% 
% lat = zeros(size(x));
% lon = zeros(size(x));
% if ~isempty(indx)
%     ce  = 2 * asin(rho(indx) / (2*radius));
%     lat(indx) = asin(D * y(indx) .* sin(ce) ./ rho(indx));
%     lon(indx) = atan2(x(indx).*sin(ce), D*rho(indx).*cos(ce));
% end