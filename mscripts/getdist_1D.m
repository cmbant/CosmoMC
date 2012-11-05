function getdist_1D(p, root, param, plotno)
if nargin < 4
    plotno=1;
end;
pts=load([root.plotroot '_p' int2str(param.n) '.dat']);
plot(pts(:,1),pts(:,2),p.lineM{plotno},'LineWidth',p.lw1);
axis([-Inf,Inf,0,1.1]);axis manual;
set(gca,'FontSize',p.axes_fontsize); 
if p.likes
ish=ishold;
hold on;
pts=load([root.plotroot '_p' int2str(param.n) '.likes']);
plot(pts(:,1),pts(:,2),p.lineL{plotno},'LineWidth',p.lw1);
xlabel(param.label,'FontSize',p.lab_fontsize);
set(gca,'ytick',[]);
if ~ish
hold off;
end
end
