

function plot_chains(root,p1,p2,nchain)
clf;
c=colormap;
for i=1:nchain
 n=load([root '_' int2str(i) '.txt']);
 plot(n(:,p1+2),n(:,p2+2),'Color',c(round(i*size(c,1)/nchain),:))
 hold on;
end
hold off;
