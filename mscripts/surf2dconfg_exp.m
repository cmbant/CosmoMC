function hand=surf2dconfg_exp(xvals,yvals,expplt,clevels,ni,zlevelstep,color);
% based on plot2dconfg_exp to do 3d surface plots
% modified from plot2dconfgeneral just uses prob instead of -logprob
% plt=prob

% do some checks
minplt=min(min(expplt));
if (max(max(expplt))==minplt)
  (['Minimum and maximum of expplt are the same! = ',num2str(minplt)])
  return; 
end;

% renormalise expplt and set threshold
expplt=expplt/max(max(expplt)); % norm is arbitrary anyway
expplt(find(expplt==0))=1e-2;

%% interpolate and find contour levels
%[pltlevels plti xlabsi ylabsi]=conf2dlevelz...
%    (xlabs,ylabs,expplt',clevels,ninterp,zlevelstep);
%pltlevels=-log(pltlevels);
%xlabsi=squeeze(xlabsi); 
%ylabsi=squeeze(ylabsi); 

% plot it
% hand=surfl(xlabsi,ylabsi,plti');

nxvals=length(xvals);
nyvals=length(yvals);

% interpolate in 2d it to make it prettier
xvalsi=interp1(1:nxvals,xvals,1:1/ni:nxvals);
yvalsi=interp1(1:nyvals,yvals,1:1/ni:nyvals);
[xmesh,ymesh]=meshgrid(xvals,yvals);
[xmeshi,ymeshi]=meshgrid(xvalsi,yvalsi);
ximi2d=interp2(xmesh,ymesh,expplt,xmeshi,ymeshi);
% make it look as pretty as possible.. use a spline to interpolate
%ximi=interp3(xmesh,ymesh,zmesh,xim,xmeshi,ymeshi,zmeshi,'spline');
% takes too long and too much memory

% use volume visualisation methods instead
zvalsi=0:0.01:1;
ximi=zeros(length(yvalsi),length(xvalsi),length(zvalsi));
for ix=1:length(xvalsi)
    for iy=1:length(yvalsi)
        % do ix and iy the bizzare way round needed by isosurf functions
        ximi(iy,ix,:)=ximi2d(ix,iy)-zvalsi;
    end
end
%plot3dconf(xvalsi,yvalsi,zvalsi,xim,0.68,ninterp,0)

% need a 3d mesh now
[xmeshi,ymeshi,zmeshi]=meshgrid(xvalsi,yvalsi,zvalsi);

levtmp=0;

% draw the objects
hpatch=patch(isosurface(xmeshi,ymeshi,zmeshi,ximi,levtmp),...
    'FaceColor','red','EdgeColor','none');
hendcaps=patch(isocaps(xmeshi,ymeshi,zmeshi,ximi,levtmp),...
      'FaceColor','interp','EdgeColor','none');
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
%hlight(2)=light('Position',[2*v(2)-v(1) 2*v(4)-v(3) 2*v(6)-v(5)],'style','local');
hlight(3)=light('Position',[v(1)-(v(2)-v(1)) v(4) (v(5)+v(6))/2],'style','local');

hand=[hpatch hendcaps hlight];

