%% pretty green colormap
cmaptmp=colormap
cmap(:,1)=1.0*((find(ones(length(cmaptmp),1))-1)/(length(cmaptmp)-1))
cmap(:,2)=1.0
cmap(:,3)=0.2+0.8*((find(ones(length(cmaptmp),1))-1)/(length(cmaptmp)-1))
colormap(cmap)




%%% pretty green colormap
%ncolsteps=64; i=ncolsteps/64; clear cmap
%cmap(:,1)=[0.9*([1:64*i]-1)/(64*i-1)]';
%cmap(:,2)=[0.9*ones(1,64*i)]';
%cmap(:,3)=[0.2+0.7*([1:64*i]-1)/(64*i-1)]';
%colormap(cmap)

%%% pretty green colormap
%ncolsteps=64; i=ncolsteps/64; clear cmap
%cmap(:,1)=[([1:64*i]-1)/(64*i-1)]';
%cmap(:,2)=[ones(1,64*i)]';
%cmap(:,3)=[0.2+0.8*([1:64*i]-1)/(64*i-1)]';
%colormap(cmap)
