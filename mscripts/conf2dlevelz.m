function [zlevels, z, x, y]=conf2dlevelz(xcoarse,ycoarse,zcoarse,clevels,ninterp,zlevelstep)
% zlevelstep=100; was assumed in conf2dlevel
% the only difference between this prog and conf2dlevel is that
% you can specify zlevelstep in this one

% ycoarse=matrx of probabilities
% xcoarse= abscissa for y (must be evenly spaced)
% clevels=[0.68,0.90,0.95] 
% I assume contour levels are increasing
% ninterp is number of points to interpolate with
% ylevelstep

meth='linear';
%meth='cubic';

zlevels=zeros(size(clevels));

[nycoarse,nxcoarse]=size(zcoarse);
nx=nxcoarse*ninterp;
ny=nycoarse*ninterp;
ystep=(ycoarse(nycoarse)-ycoarse(1))/ny;
xstep=(xcoarse(nxcoarse)-xcoarse(1))/nx;
if (ninterp>1)
  x=xcoarse(1)+(find(ones(nx,1))-1)*xstep;
  y=ycoarse(1)+(find(ones(ny,1))-1)*ystep;
%  (['Interpolating using ',meth,' with ',num2str(ninterp),' points in each direction...'])
  z=interp2(xcoarse,ycoarse,zcoarse,(ones(ny,nx)*diag(x)),diag(y)*ones(ny,nx),meth);
else
  x=xcoarse;
  y=ycoarse;
  z=zcoarse;
end

% contour(x,y,z)
% find maximum in data
zmax=max(max(z));
[jmax,imax]=find(z==max(max(z)));

% normalise likelihood function - integrate over whole range
% get rid of the NaNs
sz=size(z);
for i=1:sz(1) 
    for j=1:sz(2) 
        if (isnan(z(i,j))) 
            z(i,j)=0;
        end; 
    end; 
end
norm=sum(sum(z)) *xstep*ystep;
% step through levels in z, evaluating confidence in each 
zlevel=zmax;
conf=0;
% I assume contour levels are increasing
for i=1:length(clevels)
    fprintf(1,['Finding contour level ',num2str(i),' ...\n']);
    while (conf<clevels(i)) % so it stops when conf is just greater than clevels(i)
        zlevel=zlevel-zmax/zlevelstep;
        %% make a vector to integerate which is the same as z, but zero where z.lt.lim
        tempmat=z;
        tempmat(find(z<zlevel))=0.;
        %% integrate it, normalise, compare with confidence level reqd
        conf=(sum(sum(tempmat)) *xstep*ystep)/norm;
    end;
    while ((conf==1.0)&(zlevelstep<10000))
        if (i>=2) zlevel=zlevels(i-1); end
        if (i==1) zlevel=zmax; end
        zlevelstep=zlevelstep*2;
        (['Increasing zlevelstep to ',num2str(zlevelstep),' in conf2dlevel to improve accuracy'])
        conf=0;
        while (conf<clevels(i))
            zlevel=zlevel-zmax/zlevelstep;
            %% make a vector to integerate which is the same as z, but zero where z.lt.lim
            tempmat=z;
            tempmat(find(z<zlevel))=0.;
            %% integrate it, normalise, compare with confidence level reqd
            conf=(sum(sum(tempmat)) *xstep*ystep)/norm;
        end;
    end;
    fprintf(1,'Contour level %i is actually the %5.4g per cent contour\n',...
        i,conf*100);
    zlevels(i)=zlevel;
end;

if (any(zlevels<0)) 
    (['!!!! ERROR in conf2dlevel !!!! a zlevel is less than zero'])
end

%zlevel
%contour(x,y,z,[zlevel,zlevel])
%contour(xcoarse,ycoarse,zcoarse,[zlevel,zlevel])



