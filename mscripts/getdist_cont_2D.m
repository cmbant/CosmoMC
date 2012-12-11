function [C, h]=getdist_cont_2D(p, root, param1, param2, plotno)

if nargin < 5
    plotno=1;
end;

matname = [root.plotroot '_2D_' param1.n '_' param2.n];
if ~exist(matname,'file')
 matname = [root.plotroot '_2D_' param2.n '_' param1.n];
 trans=true;
else
 trans=false;
end

pts=load(matname);
if ~trans 
  pts=pts';
end

tmp = load ([root.plotroot '_p_' param1.n '.dat']);
x1 = tmp(:,1);
tmp = load ([root.plotroot '_p_' param2.n '.dat']);
x2 = tmp(:,1);

conts = load([matname '_cont']);

set(gca,'climmode','manual');
[C ,h] = contour(x1,x2,pts,conts,p.lineM{plotno},'linewidth',2);
set(gca,'Layer','top','FontSize',p.axes_fontsize);
