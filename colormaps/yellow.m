%% cyan colormap
cmap = colormap(gray);
cmap(:,3)=cmap(length(cmap):-1:1,1);
cmap(:,1)=1;
cmap(:,2)=1;
colormap(cmap);
