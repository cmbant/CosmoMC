function plot4D(varargin)
%Usage plot4D('file_root',param1,param2,param3,param_color)
%e.g. plot4D('WMAP',16,20,17,4)
%Note must have generated at least one 3D plot when running getdist

fname = varargin{1};

clf
d=load([varargin{1} '_single.txt']);
scatter3(d(:,varargin{2}+2),d(:,varargin{3}+2),d(:,varargin{4}+2),4,d(:,varargin{5}+2));

 fid = fopen(fullfile('plot_data',[fname '_params']));
 labs = textscan(fid, '%d %[^\n]');    
 fclose(fid);
 ix=find(labs{1}==varargin{2});
 xlabel(labs{2}(ix));
 ix=find(labs{1}==varargin{3});
 ylabel(labs{2}(ix));
 ix=find(labs{1}==varargin{4});
 zlabel(labs{2}(ix));
 ix=find(labs{1}==varargin{5});

axis tight; 
axis vis3d; %don't want auto-scaling during rotation

%Need this mess to stop colorbar moving
h=colorbar;
ps=get(h,'position')
set(h,'position',ps);
set(gca,'outerposition',[0 0 ps(1) 1]);

set(get(h,'ylabel'),'String',labs{2}(ix));
rotate3d on;
box on;
grid off;

