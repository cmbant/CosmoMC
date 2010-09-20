%% pretty blue colormap
cmap = colormap('gray');
cmap(:,1)=cmap(length(cmap):-1:1,1);
cmap(:,2)=cmap(:,1);
cmap(:,3)=1.0;
colormap(cmap);

