function contourf_col(varargin)
%Like contourf except takes sixth argument which is array of colours to use for filled contours

 [C h]= contourf(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
 
 if nargin>5 

 cols=varargin{6};
 if ischar(cols)    
    if strcmp(cols,'1')
     cols = [.15 .15 .5;0 1 1];
    elseif strcmp(cols,'2')
     cols = [1 0.2 0;1 .8 0];
    else 
    ht=text(0,0,'tmp','Color',cols(1));
    res1 = 0.75*get(ht,'Color')+0.25*[1 1 1];
    if (size(cols,2)==1)
     cols=[res1; 0.5*res1 + 0.5*[1 1 1]];
    else
     set(ht,'Color',cols(2));
     res2 = get(ht,'Color');
     cols=[res1;res2];
    end
    delete(ht);
    end
  end

 conts = varargin{4};
 pts = findobj(h, 'type','patch');
 cd=get(pts,'CData');
 if (size(cd,1)==1)
  cdats=[cd];   
 else
  cdats=sort(cell2mat(cd),'descend');
 end
 ncol=size(cdats,1);
 c(1)=cdats(1);
 ncont =1;
 for j=2:ncol
  if cdats(j) ~= c(1)
   c(2)=cdats(j);
   ncont=2;
  end
 end

for i=1:min(ncont,2)
  imh = findobj(h, 'type','patch','CData',c(i));
  set(imh,'FaceColor',cols(i,:));
  set(imh,'UserData','confid2D');
 end 

% for i=1:ncol
%  imh = findobj(h, 'type','patch','FaceVertexCData',conts(i));
%  set(imh,'FaceColor',cols(i,:),'AlphaDataMapping','none','FaceAlpha',alpha,'EdgeAlpha',alpha);
% end 
end
 

