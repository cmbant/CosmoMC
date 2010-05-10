function h=plot3dconf(xvals,yvals,zvals,xim,clev,ni,doendcaps);
% xim is propto exp(-loglike)

if (~exist('doendcaps')) doendcaps=1; end % do the endcaps by default

nxvals=length(xvals);
nyvals=length(yvals);
nzvals=length(zvals);

% interpolate it to make it prettier
xvalsi=interp1(1:nxvals,xvals,1:1/ni:nxvals);
yvalsi=interp1(1:nyvals,yvals,1:1/ni:nyvals);
zvalsi=interp1(1:nzvals,zvals,1:1/ni:nzvals);
[xmesh,ymesh,zmesh]=meshgrid(xvals,yvals,zvals);
[xmeshi,ymeshi,zmeshi]=meshgrid(xvalsi,yvalsi,zvalsi);
ximi=interp3(xmesh,ymesh,zmesh,xim,xmeshi,ymeshi,zmeshi);
% make it look as pretty as possible.. use a spline to interpolate
%ximi=interp3(xmesh,ymesh,zmesh,xim,xmeshi,ymeshi,zmeshi,'spline');
% takes too long and too much memory

% find the level corresponding to clev (eg. 0.68)
levtmp=confndlevelz(ximi,clev,200);
% assume that don't need super high accuracy for exact confidence leve
% .. using 500 instead of 100 takes a while
%levtmp=max(max(max(ximi)))*exp(-1.76), % 1sig for 3d Gaussian

% draw the objects
hpatch=patch(isosurface(xmeshi,ymeshi,zmeshi,ximi,levtmp),...
    'FaceColor','red','EdgeColor','none');
if (doendcaps)
  hendcaps=patch(isocaps(xmeshi,ymeshi,zmeshi,ximi,levtmp),...
      'FaceColor','interp','EdgeColor','none');
end
isonormals(xmeshi,ymeshi,zmeshi,ximi,hpatch);

% sort out some default views and lighting
view([-20 30]); 
axis tight; 
axis vis3d; 
grid on
lighting gouraud; % lighting gouraud is faster than phong
%alpha(0.5); % transparent 
material shiny;
box on
v=axis;
hlight(1)=light('Position',[v(2) v(3)-(v(4)-v(3)) (v(5)+v(6))/2],'style','local');
hlight(2)=light('Position',[2*v(2)-v(1) 2*v(4)-v(3) 2*v(6)-v(5)],'style','local');
hlight(3)=light('Position',[v(1)-(v(2)-v(1)) v(4) (v(5)+v(6))/2],'style','local');

h=[hpatch hendcaps hlight];