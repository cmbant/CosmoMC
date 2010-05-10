function [xim, x1vals,x2vals]=samp2im(x,count,ninbox)
% works in 2d only
% takes a list of samples x
% x1=x(1,:)
% x2=x(2,:)
% and calculates counts in cells for xim

if (~exist('ninbox')) ninbox=4; end
% if ninbox is less than 1 assume this is the resolution instead of the no
% in the box
% here resolution means width of bin as fraction of plot side

nd=2;
nsamp=length(x);

for i=1:nd
  xmin(i)=min(x(i,:));
  xmax(i)=max(x(i,:));
end

if (ninbox>1)
  nsteps=floor(nsamp^(1/nd) /ninbox);
else 
  nsteps=floor(1/ninbox);
end
  
xstep(1)=(xmax(1)-xmin(1))/nsteps;
x1vals=(xmin(1):xstep(1):(xmin(1)+(nsteps-1)*xstep(1)))+xstep(1)*0.5;
xstep(2)=(xmax(2)-xmin(2))/nsteps;
x2vals=(xmin(2):xstep(2):(xmin(2)+(nsteps-1)*xstep(2)))+xstep(2)*0.5;

xim=zeros(nsteps);
for isamp=1:nsamp
  i1=floor((x(1,isamp)-xmin(1))/xstep(1)*0.999)+1;
  i2=floor((x(2,isamp)-xmin(2))/xstep(2)*0.999)+1;
  xim(i1,i2)=xim(i1,i2)+count(isamp);
end

xim=xim/nsamp;
