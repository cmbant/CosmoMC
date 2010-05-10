%% 9 Level Color Scale Colormap with Mapping to Grayscale for Publications.
%%then interpolated
%
amap=colormap([0.1 0.1 0.15;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5;1 1 1]);
sc=20;
sz=(size(amap,1)-1)*sc+1
newmap=zeros(sz,3);
for i=1:size(amap,1);
 ii=(i-1)*sc+1;
 newmap(ii,:) = amap(i,:);
 if i<size(amap,1)
 for j=ii+1:ii+sc-1
  for c=1:3
  newmap(j,c) = amap(i,c) + (amap(i+1,c) - amap(i,c))/sc*(j-ii);
  end;
 end;
 end;
end;
colormap(newmap)


%colormap([0 0 0;0.05 0.05 0.15; 0.1 0.1 0.3; .15 .15 .5; 0.2 0.15 0.62; 0.25 0.15 0.58; .3 .15 .75; 0.4 0.16 0.5; 0.5 0.17 0.5;  ;.6 .2 .50;
%0.74 0.22 0.35;0.85 0.24 0.25; 1 .25 .15; 0.97 0.35 0.1; 0.93 0.43 0.05;  .9 .5 0; 0.9 0.65 0.05; .9 .75 .1; 0.9 0.83 0.3;  .9 .9 .5;1 1 1])

%colormap([.15 .15 .5;.3 .15 .75;.6 .2 .50;1 .25 .15;.9 .5 0;.9 .75 .1;.9 .9 .5])
%colormap([.15 .15 .5;.9 .5 0])

%colorbar