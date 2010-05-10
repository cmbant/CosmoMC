%%lighten the current colormap
cmaptmp=colormap;
cmap(:,1)=0.1+0.9*cmaptmp(:,1);
cmap(:,2)=0.1+0.9*cmaptmp(:,2);
cmap(:,3)=0.1+0.9*cmaptmp(:,3);
colormap(cmap);
