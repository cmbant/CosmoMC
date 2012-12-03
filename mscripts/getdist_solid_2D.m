function [C, h]=getdist_solid_2D(p, root, param1, param2, col,contcols)
%e.g. getdist_solid_2D(p,pl,nnu,omegak,'-k','gc')

matname = [root.plotroot '_2D_' int2str(param1.n) '_' int2str(param2.n)];
if ~exist(matname,'file')
 matname = [root.plotroot '_2D_' int2str(param2.n) '_' int2str(param1.n)];
 trans=true;
else
 trans=false;
end

pts=load(matname);
if ~trans 
  pts=pts';
end

tmp = load([root.plotroot '_p' int2str(param1.n) '.dat']);
x1 = tmp(:,1);
tmp = load([root.plotroot '_p' int2str(param2.n) '.dat']);
x2 = tmp(:,1);

conts = load([matname '_cont']);
shade=false;

if nargin>4
 if (nargin >5)
  shade=true;
  end;
else
 col='-k';
end

[C, h] = contour(x1,x2,pts,conts,':k');

ish=ishold;
hold on; axis manual; 

if shade
 contourf_col(x1,x2,pts,conts,col,contcols);
else
 [C, h]=contour(x1,x2,pts,conts,col,'LineWidth',2);
end

if ~ish
hold off;
end

if strcmp(get(get(gca,'xlabel'),'String'),'')
 xlabel(param1.label);
 if strcmp(get(get(gca,'ylabel'),'String'),'')
  ylabel(param2.label);
 end
end

