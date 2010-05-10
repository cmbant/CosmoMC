function smoothed=smoothitf1d(data,xsmw);
%

nx=length(data);

% double pad the data, use first and last entries to pad..
data2=zeros(1,nx*2);
i1=(floor(nx/2)+1);
i2=(3*floor(nx/2));
data2(i1:i2)=data;
data2(1:(i1-1))=data(1);
data2((i2+1):(nx*2))=data(nx);

% make smoothing kernel
b=zeros(nx,1);
for ix=1:nx
  b(ix)=exp(-(ix-(floor(nx/2)+1))^2/(2*xsmw^2));
end
norm=sum(sum(b));
b=b/norm;

% convolve
smoothed=conv(data2,b)';

% make the output the same length as the input..
i1=(((floor(nx/2)+1)+(floor(nx/2))));
i2=i1+nx-1;
smoothed=smoothed(i1:i2);

