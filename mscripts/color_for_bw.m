%Colormap that should print consistently in B & W
%Yellow to blue/black
cmap = colormap('gray');
sz=256;
p=1;
p2=round(sz/3);
cmap(p:p2,2)=1;
for i=p:p2
cmap(i,1)=1 - (i-1)/(sz/3);
cmap(i,3)=(i-1)/(sz/3);
end;

p=p2+1;
p2=round(2*sz/3);
cmap(p:p2,3)=1;
cmap(p:p2,1)=0;
for i=p:p2
cmap(i,2)=1 - (i-p)/(sz/3);
end;

p=p2+1;
p2=sz;
cmap(p:p2,2)=0;
cmap(p:p2,1)=0;
for i=p:p2
cmap(i,3)=1 - (i-p)/(sz/3);
end;

colormap(cmap);