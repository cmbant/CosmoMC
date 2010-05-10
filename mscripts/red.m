%% pretty blue colormap
cmaptmp=colormap;
cmap(:,1)=ones(length(cmaptmp),1);
cmap(:,2)=0.2+0.8*((find(ones(length(cmaptmp),1))-1)/(length(cmaptmp)-1));
cmap(:,3)=0.2+0.8*((find(ones(length(cmaptmp),1))-1)/(length(cmaptmp)-1));
colormap(cmap);
