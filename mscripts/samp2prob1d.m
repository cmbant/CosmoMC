function [p, x]=samp2prob1d(y,w,smw,margmeth)
% SAMP2PROB1D Probability in 1d from a list of samples
%   [P X]=SAMP2PROB1D(Y,W,S) turns a list of samples Y (with weights W)
%   into a probability distribution P defined at points X. Uses smoothing
%   length S, where S is a fraction of XMAX-XMIN
%
%   [P X]=SAMP2PROB1D(Y,W,S,MARGMETH) uses 
% MARGMETH method to marginalise. Options are 'marg' (default)
% 'mean' = mean likelihood and 'maxi' max likelihood in pixel
% If using 'mean' then W must be the likelihood rather than weight.

% SLB 

if (~exist('margmeth')) margmeth='marg'; end

xmax=max(y);
xmin=min(y);

if (xmax~=xmin)

  nstep=100;
  % just do 100 steps every time..
  if (smw<1/nstep) fprintf(1,'Increase nstep in samp2prob1d...'); end
  % NB. due to a mess in smoothitf1d then nstep has to be even

  xstep=(xmax-xmin)/nstep;
  x=(xmin+xstep/2):xstep:(xmax-xstep/2+0.0001*xstep);

if (all(margmeth=='mean'))
  n=meanhist(y,w,x);
  elseif (all(margmeth=='maxi'))
%  fprintf(1,'Haven't writen this yet..');
else
  n=whist(y,w,x);
end
  
  % bar(x,n,'hist')
  
  p=smoothitf1d(n,smw*nstep);

else 
  p=0;
  x=xmin;
end

