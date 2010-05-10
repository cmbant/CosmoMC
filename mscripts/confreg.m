function [minus, plus]=confreg(xcoarse,ycoarse,clevel,ninterp,ansplot)
% xcoarse= abscissa for ycoarse (must be evenly spaced)
% ycoarse=vector of probabilities
% clevel=0.68 
% ninterp is number of points to interpolate with
% ylevelstep
% ansplot ='n' for no plotting

% plot(xcoarse,ycoarse,'b-+')
ylevelstep=100;

%% check that y is not flat
minus=xcoarse(1);
plus=xcoarse(length(xcoarse));
if (max(ycoarse)==min(ycoarse))
  return
end
if (isnan(max(ycoarse)))
  if (isnan(min(ycoarse)))
    return
  end
end

%% interpolate
ncoarse=length(xcoarse);
n=ncoarse*ninterp;
xstep=(xcoarse(ncoarse)-xcoarse(1))/n;
x=xcoarse(1)+(find(ones(n,1))-1)*xstep;
if (ansplot(1:1)=='c')
  y=interp1(xcoarse,ycoarse,x,'cubic');
else
  y=interp1(xcoarse,ycoarse,x);
end
% hold on; plot(x,y,'r')

% find maximum in data
ymax=max(y);
imax=find(y==max(y));

% normalise likelihood function - integrate over whole range
norm=sum(y(find(isnan(y)==0)))*xstep;
if (ansplot(1:1)~='n') 
  plot(x,(y/norm),ansplot); 
  axis([x(1), x(n), 0, 1.05*ymax/norm]);
end
if ((length(ansplot)>=4)&(ansplot(1:4)=='norm'))
  if (length(ansplot)>=5)
    plot(x,(y/ymax),ansplot(5:length(ansplot))); 
    axis([x(1), x(n), 0, 1.05]);
  else
  end
end

% step through levels in y, evaluating confidence in each 
ylevel=ymax;
conf=0;
while (conf<clevel)
  ylevel=ylevel-ymax/ylevelstep;
  %% make a vector to integerate which is the same as y, but zero where y.lt.lim
  tempvect=y;
  tempvect(find(y<ylevel))=0.;
  %% integrate it, normalise, compare with confidence level reqd
  conf=(sum(tempvect) *xstep)/norm;
end;

% start at peak and step backwards along x axis until y<ylevel
[trash,i]=min((y(1:imax)-ylevel).^2);
minus=x(i);

% start at peak and step forward along x axis until y<ylevel
[trash,i]=min((y(imax:n)-ylevel).^2);
plus=x(i+imax-1);
%y(i+imax-1)
if (ansplot(1:1)~='n')
  hold on
  plot([minus,minus],[0,1.05*ymax/norm],'b:')
  plot([plus,plus],[0,1.05*ymax/norm],'b:')
  hold off
end
%if ((length(ansplot)>=4)&(ansplot(1:4)=='norm'))
%  hold on
%  if (length(ansplot)>=5)
%    plot([minus,minus],[0,1.05],ansplot(5:length(ansplot))
%    plot([plus,plus],[0,1.05],ansplot(5:length(ansplot)))
%  else
%    plot([minus,minus],[0,1.05],'b:')
%    plot([plus,plus],[0,1.05],'b:')
%  end
%  hold off
%end






