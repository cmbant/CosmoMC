function [plot_col, plot_row]=getdist_plots_1D(p, root, params, varargin)

nparam=size(params,2);

plot_row =  round(sqrt(nparam/1.4));
plot_col = (nparam +plot_row-1)/plot_row ;

for i=1:nparam 
subplot(plot_row,plot_col,i)

getdist_1D(p,root,params{i});
for j=1:nargin-3
 hold on; 
 getdist_1D(p,varargin{j},params{i},j+1);
end;
end;

%conveneience so we don't have normally have to keep track
p.last.plot_col=plot_col;
p.last.plot_row=plot_row;

end