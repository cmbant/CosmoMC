function n=whist(y,w,x)
%
% WHIST Weighted histogram.
%     N=WHIST(Y,W,X) bins the weighted elements of Y into
%     bins with centers specified by X, returning the number
%     of elements in each container.
%
%     The x bins must be equally spaced.
%
% see also HIST

% S.L.Bridle 29 May 2002

sumw=sum(w);

binw=x(2)-x(1);

for i=1:length(x)
  is=find( (y>(x(i)-binw/2)) & (y<(x(i)+binw/2)) );
%BUG!!!  n(i)=sum(y(is).*w(is)) /sumw;
  n(i)=sum(w(is)) /sumw;
end
