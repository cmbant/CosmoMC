%% colormap for b/w plots
cmaptmp=colormap;
clear cmap
cmap=cmaptmp;
cmap(:,1)=((find(ones(length(cmaptmp),1))-1)/(length(cmaptmp)-1));
cmap(:,2)=((find(ones(length(cmaptmp),1))-1)/(length(cmaptmp)-1));
cmap(:,3)=((find(ones(length(cmaptmp),1))-1)/(length(cmaptmp)-1));
colormap(cmap);

