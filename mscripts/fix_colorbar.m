%Fix for matlab 7 colorbar/label overlap bug
function fix_colorbar(hbar,ax)

set(hbar,'location','manual','ActivePositionProperty','OuterPosition');
hbarPos = get(hbar,'OuterPosition');
%set(hbar,'OuterPosition',[0 hbarPos(2) 1 hbarPos(4)])

% What are the axes positions and margin info
axPos = get(ax,'Position');
axMargin = get(ax,'TightInset');

% Calculate and set new position of axes to accomodate colorbar and margins
newAxPos = [axPos(1),hbarPos(4)+axMargin(2)+hbarPos(2),...
axPos(3), axPos(4)+axPos(2)-hbarPos(4)-axMargin(2)-axMargin(4)-hbarPos(2)];
set(ax,'ActivePositionProperty','Position','Position',newAxPos)
