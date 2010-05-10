function hand=plot2dconfg_exp(xlabs,ylabs,expplt,clevels,ninterp,zlevelstep,color)
% modified from plot2dconfgeneral just uses prob instead of -logprob
% plt=prob

% renormalise expplt and set threshold
expplt=expplt/max(max(expplt)); % norm is arbitrary anyway
expplt(find(expplt==0))=1e-6;
plt=-log(expplt);
%plt(find(isinf(plt)))=NaN;

%%%%     look at plt, the parnx,parny plane
minplt=min(min(plt));
if (max(max(plt))==minplt)
  (['Minimum and maximum of plt are the same! = ',num2str(minplt)])
  return; end;
contst=minplt;
%axis([xlabs(1) xlabs(nvals(nx)) ylabs(1) ylabs(nvals(ny))])

% are plots currently being held?
washeld=ishold;

plttmp=log(plt-min(min(plt))+1);
plttmp(find(isnan(plttmp)))=max(max(plttmp));

if (color(1)=='c')
hand=contourf(xlabs,ylabs,plttmp', ...
    ((max(max(plttmp)))-min(min(plttmp)))*(find(ones(64,1))-1)/64.+min(min(plttmp)));
%'
  shading flat
  if (washeld==0) hold on; end
elseif (color(2)=='c')
hand=contourf(xlabs,ylabs,plttmp', ...
    ((max(max(plttmp)))-min(min(plttmp)))*(find(ones(20,1))-1)/20.+min(min(plttmp)));
%'
  shading flat
  if (washeld==0) hold on; end
end

if (color(1)~='n')
  %% get contour levels
  % clevels=[0.68,0.90,0.95];
  [pltlevels plti xlabsi ylabsi]=conf2dlevelz(xlabs,ylabs,expplt',clevels,ninterp,zlevelstep);
%  fprintf(1,' %12.8g ',pltlevels);
  pltlevels=-log(pltlevels);
  xlabsi=squeeze(xlabsi); ylabsi=squeeze(ylabsi); plti=-log(plti)';
  if ((length(color)>=3)&(color(1:3))=='qbw'&length(color)>3)
  %  contour(xlabsi,ylabsi,plti',pltlevels,deblank(color(4:length(color))))
    [c,hand]=contour(xlabs,ylabs,plt',pltlevels,deblank(color(4:length(color))));
  else
  %  contour(xlabsi,ylabsi,plti',pltlevels,'k')
    [c,hand]=contour(xlabs,ylabs,plt',pltlevels,'k');
  end
%'
end

% conf2dlevelz is the same as conf2dlevel but has an extra (the last) 
% argument which says the value to use for zlevelstep
% the value used by conf2dlevel is 100.
% xcoarse=xlabs; ycoarse=ylabs; zcoarse=exp(-plt');
% ninterp=1; zlevelstep=100

if (color(1:2)~='bw')
  if (washeld==0)
    hold off
  end
end










