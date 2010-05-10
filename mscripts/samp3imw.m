function [xim, xnvals]=samp3imw(x,count,ninbox)
% based on samp2imw to do 3d
% takes a list of samples x
% x1=x(1,:)
% x2=x(2,:)
% ... to nd
% and calculates counts in cells for xim
% xnvals(1,:)=x1vals
% ... to nd

if (~exist('ninbox')) ninbox=4; end
% if ninbox is less than 1 assume this is the resolution instead of the no
% in the box
% here resolution means width of bin as fraction of plot side

tmp=size(x);
nd=tmp(1);
nsamp=tmp(2);

for i=1:nd
  xmin(i)=min(x(i,:));
  xmax(i)=max(x(i,:));
end

if (ninbox>1)
  nsteps=floor(nsamp^(1/nd) /ninbox);
else 
  nsteps=floor(1/ninbox);
end
% do the same number of steps in each direction

for id=1:nd
  xstep(id)=(xmax(id)-xmin(id))/nsteps;
  xnvals(id,:)=(xmin(id):xstep(id):(xmin(id)+(nsteps-1)*xstep(id))) ...
      +xstep(id)*0.5;
end

% why can't this bit be general for nd? Can't think how to do it..
xim=zeros(nsteps,nsteps,nsteps);
for isamp=1:nsamp
  for id=1:nd
    ii(id)=floor((x(id,isamp)-xmin(id))/xstep(id)*0.999)+1;
  end
  % make the array the bizzare way around required by isosurface etc
  xim(ii(2),ii(1),ii(3))=xim(ii(2),ii(1),ii(3))+count(isamp);
end

%xim=xim/max(max(max(xim))); % so peak is at 1
