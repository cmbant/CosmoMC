%% purple colormap
cmap = colormap(gray);
cmap(:,2)=cmap(length(cmap):-1:1,1);
cmap(:,1)=1;
cmap(:,3)=1;
colormap(cmap);
