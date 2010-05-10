%% invert the colormap
cmaptmp=colormap;
cmap(1:length(cmaptmp),1)=cmaptmp(length(cmaptmp):-1:1,1);
cmap(1:length(cmaptmp),2)=cmaptmp(length(cmaptmp):-1:1,2);
cmap(1:length(cmaptmp),3)=cmaptmp(length(cmaptmp):-1:1,3);
colormap(cmap);
