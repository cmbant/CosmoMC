%% green colormap
cmap = colormap(gray);
cmap(:,1)=cmap(length(cmap):-1:1,1);
cmap(:,2)=1;
cmap(:,3)=cmap(:,1);
colormap(cmap);
