function showContoursBehind(varargin)

ax = findobj(gcf,'Type','axes');

for j=1:size(ax,1);

t = findobj(ax(j),'Type','hggroup');

for i=2:size(t,1)

newc=copyobj(t(i),ax(j));

p = findobj(newc,'Type','patch');
if nargin>0
 set(p,'LineStyle',varargin{1});
else
 set(p,'LineStyle','--');
end
set(p,'FaceColor','none');

end
end