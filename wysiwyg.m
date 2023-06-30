function wysiwyg
% WYSIWYG -- Resize figure on screen according to papersize.
%      
%       This function is called with no args and merely
%       changes the size of the figure on the screen to equal
%       the size of the figure that would be printed, 
%       according to the papersize attribute.  Use this function
%       to give a more accurate picture of what will be 
%       printed.

%       Dan(K) Braithwaite, Dept. of Hydrology U.of.A  11/93

pause(0.01) %Seemed to be issues if called too fast after changing some things AGR,10/11/2013
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
pos(3:4) = ppos(3:4);
% pos(1:2) = [1 -1];
pos(1:2) = [0.1 0.6];
% pos(1:2) = [1 6];
set(gcf,'position',pos);
set(gcf,'units',unis);
  
