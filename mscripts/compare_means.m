%Makes horizontal bar chart of the difference between means 
%(from .margestats files) in units of the errors bars
%Useful for accuracy testing/seeing changes between analysis options.
%Results are relative to the first input e.g.
% compare_means('high_accuracy','test','test2');
%
function compare_means(varargin)

%varargin={'post2','post1','post12'};
%nargin=size(varargin,2);

ext='.margestats';
data={};
errs=[];
names=[];
manes=[];
for i=1:nargin

    ncols=0;
    nvar=0;    
    fid=fopen([varargin{i} ext]);
    header=fgetl(fid);
    while true 
        tmp=[];
        [tmp, ncols] = fscanf(fid, '%e');
        if (ncols<1)
             break;
        end;   
        nvar=nvar+1;
        means(i,nvar)=tmp(2);
        if (i==1)
       %  names{nvar}=['$' fgetl(fid) '$'];
         names{nvar}=[fgetl(fid)];

         errs(nvar)=tmp(3);
        else
           fgetl(fid); 
        end;
    end;
    fclose(fid);

end;

compare=means(1,:);

fracs=[];
for i=2:nargin
    fracs(:,i-1)= (means(i,:)-compare)./errs;
end;

barh(fracs);
for i=1:nvar
text(0.05,(i)/(nvar+1),names{i},'units','normalized','fontsize',14);
end;
set(gca,'ytick',[]);
title('difference in posterior mean / error bar');
legend(varargin{2:nargin},'location','Best');

%set(gcf, 'PaperUnits','inches');
%set(gcf, 'PaperPosition',[ 0 0 6.8 9]);
%print -depsc2 'bar_errors.eps';