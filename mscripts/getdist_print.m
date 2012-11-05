function getdist_print(p, fname, wid, high)

if nargin<3
 if isfield(p.last,'plot_col')
    wid = p.last.plot_col;
    high=p.last.plot_row;
 else
    wid=1;
    high=1;
 end;
end

set(gcf, 'PaperUnits','inches');
x=wid*p.plot_size_inch; y=high*p.plot_size_inch;
set(gcf, 'PaperPosition',[0 0 x y]); set(gcf, 'PaperSize',[x y]);
print(p.print,[fname '.' p.printex]);