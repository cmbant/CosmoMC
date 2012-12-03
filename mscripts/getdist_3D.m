function p=getdist_3D(p, root, param1, param2, param3, varargin)

if (~isfield(root.last,'single'))
root.last.single = load([root.distroot '_single.txt']);
end
colormap(p.colormap);
h=scatter(root.last.single(:,param1.n+2),root.last.single(:,param2.n+2),3,root.last.single(:,param3.n+2));
xlabel(param1.label,'FontSize',p.lab_fontsize);
ylabel(param2.label,'FontSize',p.lab_fontsize);
set(gca,'FontSize',p.axes_fontsize); ax = gca;
hbar = colorbar('horiz');axes(hbar);
xlabel(param3.label,'FontSize',p.lab_fontsize);
set(gca,'FontSize',p.axes_fontsize);
axes(ax);

if nargin>5
   %compare 
   ish=  ishold;
   hold on;
   for i=1:nargin-5
   getdist_cont_2D(p,varargin{i},param1,param2,i); 
   end
   if ~ish
       hold off
   end
end
