function confid2D(varargin)
%Usage confid2D('file_root',param1,param2,line_color,contour_colors)
%e.g. confid2D('lya',8,17,'-k','gc')

fname = varargin{1};
p1=int2str(varargin{2});
p2=int2str(varargin{3});

matname = [fname '_2D_' p1 '_' p2];
if ~exist(fullfile('plot_data',matname),'file')
 matname = [fname '_2D_' p2 '_' p1];
 trans=true;
else
 trans=false;
end

pts=load (fullfile('plot_data',matname));
if ~trans 
  pts=pts';
end

tmp = load (fullfile('plot_data',[fname '_p' p1 '.dat']));
x1 = tmp(:,1);
tmp = load (fullfile('plot_data',[fname '_p' p2 '.dat']));
x2 = tmp(:,1);

if false
 pts=interp2(pts,2,'cubic');
 x1=interp1(x1,1:0.25:size(x1))';
 x2=interp1(x2,1:0.25:size(x2))';
end

conts = load (fullfile('plot_data',[matname '_cont']));

contcols=[];
shade=false;

if nargin>3
 col= varargin{4};
 if (nargin >4)
  shade=true;
  contcols=varargin{5};
  end;
else
 col='-k';
end

[C h] = contour(x1,x2,pts,conts,':k');
%set(h,'LineWidth',lw1);
ish=ishold;
hold on; axis manual; 

if shade
 contourf_col(x1,x2,pts,conts,col,contcols);
else
 [C h]=contour(x1,x2,pts,conts,col,'LineWidth',2);
% [C h]=contour(x1,x2,pts,conts,col,'LineWidth',2);
end

if ~ish
hold off;
end

if strcmp(get(get(gca,'xlabel'),'String'),'')
 fid = fopen(fullfile('plot_data',[fname '_params']));
 labs = textscan(fid, '%d %[^\n]');    
 fclose(fid);
 ix=find(labs{1}==varargin{2});
 xlabel(labs{2}(ix));
 if strcmp(get(get(gca,'ylabel'),'String'),'')
  ix=find(labs{1}==varargin{3});
  ylabel(labs{2}(ix));
 end
end

%set(gcf, 'PaperUnits','inches');
%set(gcf, 'PaperPosition',[ 0 0 9 8]);
%print -depsc2 'Confid2D.eps';
