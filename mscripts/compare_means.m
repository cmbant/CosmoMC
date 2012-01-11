%Makes horizontal bar chart of the difference between means 
%(from .margestats files) in units of the errors bars
%Useful for accuracy testing/seeing changes between analysis options.
%Results are relative to the first input e.g.
% compare_means('high_accuracy','test','test2');
% Also optional parameters for legend captions, to use latex, or compare errors e.g.
% compare_means('highacc','default','lowacc','legend',{'High accuracy','Default','Low accuracy'},'latex');

function compare_means(varargin)
%varargin={'post2','post1','post12'};
%nargin=size(varargin,2);


ext='.margestats';
data={};
errs=[];
names=[];
manes=[];
nargs=nargin;
interpreter='tex';
compare_errors=false; %compare means by default
for i=1:nargin
    if (strcmp(varargin{i},'legend'))
      legend_labels=  varargin{i+1}; 
      nargs=min(i-1,nargs);
    elseif (strcmp(varargin{i},'latex'))
        interpreter='latex';    
        nargs=min(i-1,nargs);
    elseif (strcmp(varargin{i},'errors'))    
        compare_errors=true;
        nargs=min(i-1,nargs);
    end;    
end;
if ~exist('legend_labels')
 legend_labels={varargin{2:nargs}};
end; 
for i=1:nargs
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
        if (compare_errors) 
         means(i,nvar)=tmp(3);
        else
         means(i,nvar)=tmp(2);
        end;
        if (i==1)
         name=[fgetl(fid)];
         if (strcmp(interpreter,'latex'))  
          if strcmp(name(1:min(3,length(name))),'log')
             name=['$\' name '$'];
          else
             name=['$' name '$'];
          end          
         end;
         names{nvar}=name;
         errs(nvar)=tmp(3);
        else
           fgetl(fid); 
        end;
    end;
    fclose(fid);

end;

compare=means(1,:);

fracs=[];
for i=2:nargs
    fracs(:,i-1)= (means(i,:)-compare)./errs;
end;

barh(fracs);
left=min(min(fracs));
right=max(max(fracs));
xlim([left-0.06 right+0.05]);
for i=1:nvar
text(left-0.05,i,names{i},'fontsize',14,'interpreter',interpreter);
end;
ticks=get(gca,'xtick');
if (ticks(2)~=0) 
set(gca,'xtick',ticks(2:size(ticks,2)));
end;
set(gca,'XTickMode','manual');
set(gca,'ytick',[]);
if (compare_errors)
title('fractional difference in error bar','fontsize',13);
else
title('difference in posterior mean / error bar','fontsize',13);
end;
legend(legend_labels,'location','Best');

print('-depsc2', [varargin{1} '_bar_errors.eps']);
