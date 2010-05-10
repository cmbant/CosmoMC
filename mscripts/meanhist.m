function n=meanhist(y,w,x)
%
% MEANHIST Mean likelihood in each cell
%     N=MEANHIST(Y,W,X) bins the weighted elements of Y into
%     bins with centers specified by X, returning the mean of
%     W in each container. (eg. w is likelihood)
%
%     The x bins must be equally spaced.
%
% see also HIST

% S.L.Bridle 29 May 2002

sumw=sum(w);

binw=x(2)-x(1);

for i=1:length(x)
  is=find( (y>(x(i)-binw/2)) & (y<(x(i)+binw/2)) );
  if (length(is)>0)
    n(i)=mean(exp(-(w(is)-min(w)))); % need to normalise so don't get zero..
  else
    n(i)=0;
  end
end
% don't think we are interested in the normalisation?
n=n/max(n);